import re     # Generate cartesian products (w-mers)
import random # Generate random sequence q
import itertools


def find(substring, string):
    '''
    Return starting indices of all occurences
    of substring in string with overlap
    '''
    return tuple(
        match.start()
        for match in re.finditer(f'(?={substring})', string)
    )


def random_str(charset, length):
    '''
    Generate random string from charset of fixed length
    '''
    return ''.join(random.choice(charset) for _ in range(length))


def chunks(iterable, size=10):
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))