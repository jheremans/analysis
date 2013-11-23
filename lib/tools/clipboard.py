from Tkinter import Tk
import numpy as np

def copy_str(str):
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(str)
    r.destroy()

def copy_array(arr):
    to_str=''
    for a in arr:
        to_str=to_str+str(a)+'\t'
    copy_str(to_str)

def copy_array_2d(arr):
    to_str=''
    for r in arr:
        for c in r:
            to_str=to_str+str(c)+'\t'
        to_str=to_str+'\n'
    copy_str(to_str)