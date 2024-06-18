""" Some simple search/replace functions. """
import shutil
import tempfile
import re
import io
import os


def open_if_needed(file, mode="r"):
    """ Open a file if the passed object is not already an open file. """
    if isinstance(file, io.TextIOBase):
        if file.mode == mode:
            return file
        file.close()
        file = open(file.name, mode)
    return open(file, "r")


def go_to_keyword(file, regex):
    """ Read through a file unil a line matches the given regular expression. """
    cfile = open_if_needed(file)
    reg = re.compile(regex)
    for line in cfile:
        if reg.search(line):
            return cfile, line
    raise ValueError('Requested keyword not found in file: ' + regex)


def check_in_file(file, regex):
    """ Check if a regex is present in a file. """
    try:
        cfile = go_to_keyword(file, regex)[0]
        cfile.close()
        return True
    except ValueError:
        return False


def search_file(file, search, after=0, max_res=None, close=True, stop_at=None, matching_only=False):
    """ Find all occurences of a regex in file. """
    cfile = open_if_needed(file)
    search_reg = re.compile(search)
    if stop_at is not None:
        stop_reg = re.compile(stop_at)
    values = []
    n_hit = 0
    for line in cfile:
        if stop_at is not None:
            if stop_reg.search(line):
                break
        line_res = search_reg.search(line)
        if line_res:
            n_hit = n_hit + 1
            if matching_only:
                values.append(line_res.groups() or line_res.group())
            elif after == 0:
                values.append(line.rstrip())
            else:
                for _ in range(after):
                    line = next(cfile).rstrip()
                    values.append(line.rstrip())
            if max_res is not None:
                if n_hit >= max_res:
                    break
    if close:
        cfile.close()
    if not values:
        raise ValueError(f"No matches for {search_reg.pattern} in file {cfile.name}")
    return values


def split_columns(in_list, col=None, split=str.split, convert=lambda x:x):
    """ Split all sublists of split_list and take specified columns. """
    if col is None:
        split_list = [[convert(sval) for sval in split(val)] for val in in_list]
    elif isinstance(col, int):
        split_list = [convert(split(val)[col]) for val in in_list]
    elif isinstance(col, list):
        split_list = [[convert(split(val)[ci]) for ci in col]for val in in_list]
    return split_list


def replace_inplace(file_name, search, replace):
    """ Do a simple in-place regex search/replace. """
    reg = re.compile(search)
    tot_n_subs = 0
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
        with open(file_name, "r") as cfile:
            for line in cfile:
                line, n_subs = reg.subn(replace, line)
                tot_n_subs += n_subs
                tmp_file.write(line)
        tname = tmp_file.name
    shutil.move(tname, file_name)
    return tot_n_subs


def replace_cols_inplace(file_name,
                         array,
                         keyword,
                         cols=None,
                         split=str.split,
                         max_replace=1,
                         skip=0,
                         val_format=" {:20.12f} "):
    """ Replace columns after a given keyword with values from array. """
    reg = re.compile(keyword)
    n_found = 0
    if cols is None:
        cols = list(range(len(array[0])))
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
        with open(file_name, "r") as cfile:
            for line in cfile:
                if not reg.search(line) or n_found >= max_replace:
                    tmp_file.write(line)
                else:
                    n_found = n_found + 1
                    tmp_file.write(line)
                    for i in range(skip):
                        line = next(cfile)
                        tmp_file.write(line)
                    for vals in array:
                        line = split(next(cfile).rstrip())
                        for i, col in enumerate(cols):
                            line[col] = val_format.format(vals[i])
                        line = ''.join(line)
                        tmp_file.write(line + '\n')
        tname = tmp_file.name
    shutil.move(tname, file_name)


def remove(path):
    """ Remove a file/directory if it exists. """
    try:
        shutil.rmtree(path)
    except NotADirectoryError:
        os.remove(path)
    except FileNotFoundError:
        pass


def fortran_double(arg):
    """ Convert fortran double format to float."""
    if isinstance(arg, bytes):
        arg = arg.decode()
    if isinstance(arg, str):
        arg = arg.replace("D", "E")
    return float(arg)
