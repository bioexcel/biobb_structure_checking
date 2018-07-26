#! /usr/bin/python3

"""
 structureCheck: Command line replacement for MDWeb structure check    
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

import sys
import os


def main():
    print ("Hellp")
    if len(sys.argv) == 0:
        help('general')
        sys.exit()
    if sys.argv[0] == '--help':
        help('general')
        sys.exit()
        
    
if __name__ == "__main__":
    main()
