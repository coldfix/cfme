import sys
import numpy.random


def main(num_samples, dimension):
    fmt = "{:8.5f}".format
    x = numpy.random.normal(size=(num_samples, dimension))
    x /= numpy.linalg.norm(x, axis=1)[:, numpy.newaxis]
    for v in x:
        print("[ {} ]".format(" ".join(fmt(c) for c in v)))


if __name__ == '__main__':
    main(*map(int, sys.argv[1:]))
