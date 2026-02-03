# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 14:15:48 2025

@author: tomke
"""

import struct as st
import decimal as dc


def exactly(num):
    exact_value = dc.Decimal(num)
    return exact_value

def num2float32(num):
    #!f indicates 32-bit
    binary32 = st.pack('!f', num)
    bits = ''.join(f'{byte:08b}' for byte in binary32)
    position = [1,9]
    for ind in sorted(position, reverse=True):
        bits = bits[:ind] + ' ' + bits[ind:]
    return bits
    
def num2float64(num):
    #!d indicates 64-bit
    binary64 = st.pack('!d', num)
    bits = ''.join(f'{byte:08b}' for byte in binary64)
    position = [1,12]
    for ind in sorted(position, reverse=True):
        bits = bits[:ind] + ' ' + bits[ind:]
    return bits
   
if __name__ == '__main__':
    x = 0.15625
    print(num2float32(x))
    print(num2float64(x))


