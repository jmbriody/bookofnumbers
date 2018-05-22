from collections import defaultdict
import math
import asyncio
# from profilehooks import coverage
import timeit

SIEVE_LIMIT = 1000000

'''
9007199254740881  (the largest prime under 2^53)
9007199254740997  (the smallest prime over 2^53)
9007199254740991  = 6361*69431*20394401
9007199254740992  = 2^53
9007199254740993  = 3*107*28059810762433
9999986200004761     = 99999931*99999931
99999989237606677    = 316227731*316227767
999999866000004473   = 999999929*999999937
9999999942014077477  = 3162277633*3162277669
99999980360000964323 = 9999999017*9999999019
99999980380000962361 = 9999999019*9999999019
99994452415200570527 = 4641503*4641503*4641503
99999999999999999999 = 3*3*11*41*101*271*3541*9091*27961
10000000000000000001 = 11*909090909090909091
10000000000000099
10,000,000,000,000,099
178571428571429
10000000000000024

'''
# Start with basic sieve variant. You can find faster and more efficient ones at places
# like stackoverflow--but I went with simple and straightforward. Ultimate goal
# is not dealing with huge primes but finding smaller primes and getting factors
# for smallish numbers.
def primes_sieve(limit):
    '''Prime Sieve (variant). Returns a list of primes upto limit.

    Rather than marking non-primes in a list it adds non-primes to a defaultdict
    and skips testing items already added. Once that is completed it subtracts a
    set of non-prime keys from a set of all numbers to limit and returns a sorted
    list.
    '''
    limitplus = limit+1
    notprimes = defaultdict(int)

    # add non-primes to notprimes defaultdict
    for i in range(2, limitplus):
        if i not in notprimes:
            temp_dict = {n: True for n in range(i + i, limitplus, i)}
            notprimes.update(temp_dict)

    # Set of notprimes keys and allnumbers upto limit
    notprime = set(notprimes.keys())
    allnumbers = set(range(2, limitplus))

    return sorted(list(allnumbers - notprime))

# Create global PRIME_LIST. On my relatively wimpy laptop "SIEVE_LIMIT" of
# 5,000,000 (348,513) primes is moderately OKish speedwise (5 to 10 seconds).
# Larger lists will of course take more memory.
# For most cases a list upto the sqrt of the largest prime factor should be good
PRIME_LIST = primes_sieve(SIEVE_LIMIT)
NEW_PRIMES = defaultdict(int)

def reset_prime_list(newlength):
    global SIEVE_LIMIT
    SIEVE_LIMIT = newlength
    global PRIME_LIST
    PRIME_LIST = primes_sieve(SIEVE_LIMIT)

def _is_prime(k):
    end = int((k**0.5) + 1)
    for j in PRIME_LIST:
        if j > end:
            return True
        elif k % j == 0:
            return False


# Probably need to break this up into function calls.
#
# @coverage
def prime_factors(n):
    original = n
    end = int((n**0.5) + 1)
    result = []

    if int(n) in PRIME_LIST or int(n) in NEW_PRIMES: # 97
        result.append((n, 1))
        n = 1


    for x in PRIME_LIST: # 96
        if x > end or n <= 1:
            break
        else:
            count = 0
            while n % x == 0:
                count += 1
                n /= x
            if count > 0:
                result.append((x, count))

    if n <= 1:
        return result

    # Check if remainder from above is prime
    if n in PRIME_LIST: # 100023
        result.append((int(n), 1, "PL"))
        n /= n

    if n > 1 and _is_prime(n): # 100020
        NEW_PRIMES[n] += 100
        result.append((int(n), 1, "B"))
        n = 1

    if n > 1: # and n not in PRIME_LIST and n not in NEW_PRIMES: # and n > end:
        current = PRIME_LIST[-1] + 2
        while n > 1 and current < (n**0.5 + 1):
            count = 0

            while n % current == 0:
                count += 1
                n //= current
            if count > 0:
                print("D", current, count, n)
                NEW_PRIMES[current] += 1
                result.append((current, count, "D"))

            current += 2
        if n > 1 and n not in NEW_PRIMES:
            NEW_PRIMES[n] += 10
            result.append((int(n), 1, "C"))
            print("Current", current)
            n = 1
        elif n in NEW_PRIMES:
            result.append((int(n), 1, "R"))
    elif n > 1 and n != original:
        result.append((int(n), 1, "E"))

    return result

# @coverage
def factors(n):
    primef = prime_factors(n)
    tempresult = [[1]]
    resset = set()

    if len(primef) == 0:
        return [1]

    # First take all prime factor tuples and generate prime ** [0 through instances]
    for x in range(0, len(primef)):
        tfactor = primef[x][0]
        expanded = [tfactor**q for q in range(1, primef[x][1] + 1) ]
        tempresult.append(expanded)
    # print(tempresult)

    resset = set(tempresult[0])
    
    # I can never find my jelly pens with twin 2.5 year-olds in the house
    # Need to take first list from factor list and convert to set
    # Iteratively multiply following lists from factor list and make them
    #    the new set
    tresult = []
    for x in range(1, len(tempresult)):
        for item in resset:
            tresult.extend([item * r for r in tempresult[x]])
        #print("resset:", resset)
        resset |= set(tresult)
        #print("tresult:", tresult)
        #print("temp", x, tempresult[x])

    return sorted(list(resset))[0:-1]   # Drop last item--the number itself

def number_type(n):
    nfactors = factors(n)
    sum_factors = sum(nfactors)
    difference = sum_factors - n
    if len(nfactors) == 1:
        return ("PRIME", difference)
    elif difference > 0:
        return ("ABUNDANT", difference)
    elif difference == 0:
        return ("PERFECT", difference)
    else:
        return ("DEFICIENT", difference)

def fibonacci(n):
    t = 1.618
    s5 = math.sqrt(5)
    fresult = ((t**n) - ((-1/t)**n)) / s5
    return math.ceil(fresult)

def factors2(n):
    return set(x for tup in ([i, n//i] 
                for i in range(1, int(n**0.5)+1) if n % i == 0) for x in tup)

if __name__ == "__main__":
    print(len(PRIME_LIST))
#    startt = timeit.default_timer()
#    for x in range(10000000, 10004000):
#        factors(x)
#    print(timeit.default_timer() - startt)

#    startt = timeit.default_timer()
#    for x in range(10000000, 10004000):
#        factors2(x)
#    print(timeit.default_timer() - startt)

    
#    reset_prime_list(1000)
#    print(len(PRIME_LIST))
#    prime_factors(100023)
#    for x in range(100000010, 100000030):
#        print(x, prime_factors(x))
#    print(NEW_PRIMES)
#    print(prime_factors(100000027))
#    print(prime_factors(66413))
