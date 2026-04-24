#!/bin/sh
'''exec' "path_to_/scripts/sgb_to_gtdb_profile.py" "$@"
' '''
# -*- coding: utf-8 -*-
import re
import sys
from metaphlan.utils.sgb_to_gtdb_profile import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
