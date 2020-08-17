#!/usr/bin/env python3

import re
import sys
import os
import os.path
import argparse
import base64
import hashlib
import pickle

from collections import defaultdict

# This matches all require statements:
# 1st group is the condition
# if this is the first time errorcodes.py runs:
#   2nd group is the message
#   3rd, 4th groups are None
# otherwise:
#   2nd group is the errorcode from prev run
#   3rd group is the whole !ERROR comment
#   4th group is the message
PATTERN = re.compile(
    r'require\s*\(\s*(.+?)\s*,\s*"(.*?)"\s*\)\s*;( /\* !ERROR: (.*?) \*/)?',
    flags=re.DOTALL)

TEMPLATE = 'require({}, "{}"); /* !ERROR: {} */'
GAS_PER_BYTE = 200

# Colors for pretty logging
CYAN = '\x1b[36m'
YELLOW = '\x1b[33m'
MAGENTA = '\x1b[35m'
RESET = '\x1b[0m'


def prefix(a, b):
    '''
    Utility function that returns shortest prefix
    of string `a` that is not in string `b`
    '''
    pref = []
    for i, j in zip(a, b):
        pref.append(i)
        if i != j:
            return ''.join(pref)


class ErrorCoder:
    def __init__(self):
        self.hashes = {}
        self.saved = 0
        self.found = False
        # the following are assigned by the argument parser
        self.recursive = False
        self.quiet = False
        self.white = False
        self.check = False
        self.code = None
        self.paths = []

    @staticmethod
    def hash(key):
        ''' Return base32(sha1(pickle(key))) '''
        serialized = pickle.dumps(key)
        hashed = hashlib.sha1(serialized).digest()
        encoded = base64.b32encode(hashed).decode()
        return encoded

    def log(self, *args, **kwargs):
        if self.quiet:
            return
        if self.white:  # replace all colors with empty strings
            args = map(lambda x: re.sub('\x1b\\[\\d\\d?m', '', x), args)
        print(*args, **kwargs)

    def process_paths(self, paths, callback):
        ''' Apply `callback` to all .sol files in `paths` '''
        for path in paths:
            if not os.path.exists(path):
                raise OSError(f"Invalid path {path}")
            if os.path.isfile(path) and \
               os.path.splitext(path)[1] == '.sol':
                callback(path)
            elif self.recursive and os.path.isdir(path):
                new_paths = (os.path.join(path, file)
                             for file in os.listdir(path))
                self.process_paths(new_paths, callback)

    def crop_hashes(self):
        ''' Leave only shortest unique prefix of every hash '''
        inverted = {v: k for k, v in self.hashes.items()}
        hashes = sorted(self.hashes.values())
        for index, string in enumerate(hashes):
            prv = nxt = string[0]
            if index > 0:
                prv = prefix(string, hashes[index-1])
            if index < len(hashes) - 1:
                nxt = prefix(string, hashes[index+1])
            longest = prv if len(prv) > len(nxt) else nxt
            self.hashes[inverted[string]] = longest

    def scan_file(self, filename, check=False):
        '''
        Scans file for require statements and populates
        self.hashes with message hashes for further processing
        If check is True, instead checks anything should be updated
        '''
        with open(filename) as file:
            source = file.read()
        counter = defaultdict(int)
        filename = os.path.split(filename)[1]
        for match in PATTERN.finditer(source):
            message = match[4] or match[2]
            key = filename, message, counter[message]
            counter[message] += 1
            if check:
                assert self.hashes.get(key) is not None
                if self.hashes[key] != match[2]:
                    self.found = True
                    self.log("require messages in",
                             f"{YELLOW}[{filename}]{RESET}",
                             "should be updated")
                    break
            else:
                self.hashes[key] = ErrorCoder.hash(key)

    def replacer(self, filename):
        '''
        Returns a function that is called each time
        we find a require statement and want to replace it;
        Returned function returns a string that should
        replace this statement
        '''
        counter = defaultdict(int)
        filename = os.path.split(filename)[1]

        def on_match(match):
            message = match[4] or match[2]
            key = filename, message, counter[message]
            counter[message] += 1
            code = self.hashes[key]
            self.saved += len(match[2]) - len(code)
            self.log(f"{CYAN}[{filename}]",
                     f"{MAGENTA}{code}{RESET} ->",
                     f"{YELLOW}{message}{RESET}")
            return TEMPLATE.format(match[1], code, message)

        return on_match

    def rewrite_file(self, filename):
        ''' Replaces all require messages in file with their errorcodes '''
        with open(filename) as file:
            source = file.read()
        rewritten = PATTERN.sub(self.replacer(filename), source)
        with open(filename, 'w') as file:
            file.write(rewritten)

    def search_file(self, filename):
        '''
        Searches for `self.code` in `filename`
        If found, prints filename and line
        '''
        with open(filename) as file:
            lines = file.readlines()
        for index, line in enumerate(lines, 1):
            for match in PATTERN.finditer(line):
                if match[2] == self.code:
                    self.found = True
                    self.log(f"{CYAN}[{filename}:{index}]",
                             f"{YELLOW}{match[0]}{RESET}")

    def run(self):
        '''
        Does the main job
        If self.code is specified, runs in search mode
        If self.check is True, runs in check mode
        Otherwise, runs in normal mode
        '''
        if self.code:  # search for self.code in paths
            self.process_paths(self.paths, self.search_file)
            if not self.found:
                self.log(f"{MAGENTA}No matches found{RESET}", file=sys.stderr)
        else:
            self.process_paths(self.paths, self.scan_file)
            self.crop_hashes()
            if self.check:  # check if something needs to be replaced
                check_file = lambda file: self.scan_file(file, check=True)
                self.process_paths(self.paths, check_file)
                if self.found:
                    self.log(f"{MAGENTA}Run errorcodes.py to fix{RESET}",
                             file=sys.stderr)
                    exit(1)
                else:
                    self.log(f"{CYAN}All error codes are up to date!{RESET}",
                             file=sys.stderr)
            else:  # replace all messages with error codes
                self.process_paths(self.paths, self.rewrite_file)
                if self.saved != 0:
                    self.log("We saved you:",
                             f"  {self.saved} bytes",
                             f"  {self.saved * GAS_PER_BYTE} gas",
                             file=sys.stderr, sep='\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--recursive', action='store_true',
                        help='recurse into directories')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="don't print anything")
    parser.add_argument('-w', '--white', action='store_true',
                        help="don't use fancy colors for logging")
    parser.add_argument('paths', metavar='PATH', nargs='+',
                        help='paths to search for .sol files in')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('-c', '--check', action='store_true',
                       help="check mode: check for messages that should be replaced")
    modes.add_argument('-f', '--find', metavar='CODE', dest='code',
                       help="search mode: find message corresponding to CODE")
    error_coder = ErrorCoder()
    parser.parse_args(namespace=error_coder)
    error_coder.run()
