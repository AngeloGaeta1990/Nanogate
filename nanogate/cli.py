import sys
import argh
from nanogate.raw_nanogate import raw_nanogate
from nanogate.filtered_nanogate import filtered_nanogate


def main():
    """
    parsing arguments
    """
    argh.dispatch_commands([raw_nanogate, filtered_nanogate
                            ])


if __name__ == "__main__":
    sys.exit(main())
