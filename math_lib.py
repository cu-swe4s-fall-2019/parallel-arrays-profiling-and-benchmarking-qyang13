import math


def list_mean(L):
    '''
    Take a numerical array as input and return the mean of the values
    ---
    Input: L
    A numerical array
    ---
    Output: m
    A float mean of the numerical array
    '''
    if L is None:
        raise ValueError("Input list cannot be None")
    if len(L) == 0:
        raise ValueError("Input list cannot be empty")

    nums = []

    for i in L:
        if type(i) is int or type(i) is float:
            nums.append(i)
        else:
            raise ValueError("Invalid value in the input array")
        continue

    m = sum(nums)/len(nums)

    return m


def list_stdev(L):
    '''
    Take a numerical array as input and return the standard deviation
    ---
    Input: L
    A numerical array
    ---
    Output: stdv
    A float standard deviation of the numerical array
    '''
    if L is None:
        raise ValueError("Input list cannot be None")
    if len(L) == 0:
        raise ValueError("Input list cannot be empty")

    nums = []

    for i in L:
        if type(i) is int or type(i) is float:
            nums.append(i)
        else:
            raise ValueError("Invalid value in the input array")
        continue

    mean = sum(nums)/len(nums)
    stdv = math.sqrt(sum([(mean-x)**2 for x in nums]) / len(nums))

    return stdv
