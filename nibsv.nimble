# Package

version       = "0.0.2"
author        = "Zev Kronenberg"
author        = "Christopher Dunn"
author        = "Brent Pedersen"
description   = "Structural Variant nibbles"
license       = "MIT"
bin           = @["nibsv"]


# Dependencies

requires "nim >= 1.2.0", "hts", "kmer", "argparse == 0.10.1", "https://github.com/brentp/nim-kmer >= 0.2.6"
