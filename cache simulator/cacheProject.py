import sys, math, time

if (len(sys.argv) < 4):
    print("Missing argument. Correct usage: sh run_sim.sh [traceFile] [cacheSizeBytes] [cacheLineSizeBytes] [numberOfWays]")
    sys.exit()

cacheSize = int(sys.argv[2])
cacheLineSize = int(sys.argv[3])
ways = int(sys.argv[4])
filename=sys.argv[1]


cache = []
cacheLines = 0
numSets = 0
linesPerSet = 0
hitCount = 0
missCount = 0


def initCache(cacheSize, cacheLineSize, ways):
    '''Function to create cache'''
    lineCount = 0
    setNr = 0
    wayNr = 0

    global cacheLines
    global numSets
    global linesPerSet

    cacheLines = cacheSize // cacheLineSize
    numSets = cacheLines // ways
    linesPerSet = cacheLines // numSets

    for i in range(cacheLines):
        cache.append({"setNr": 0, "wayNr": 0, "tag": 0, "dataExists" : 0, "writeFlag": 0, "latestTime": -1 })

    for i in range(cacheLines):
        if (lineCount == linesPerSet):
            lineCount = 0
            setNr = setNr + 1
            wayNr = 0
        cache[i]["setNr"] = setNr
        cache[i]["wayNr"] = wayNr
        cache[i]["tag"] = -1
        cache[i]["writeFlag"] = -1
        # cache[i]["dataExists"] = 0
        # cache[i]["latestTime"] = -1

        lineCount = lineCount + 1
        wayNr = wayNr + 1



initCache(cacheSize, cacheLineSize, ways)

file = open(filename, "r")
while True:
    readLine = file.readline()
    if readLine == "" or readLine == "#eof":
        break

    accessFields = readLine.split()
    if (len(accessFields) != 3):
        continue

    hitCount = hitCount + 1

    # convert the last elemtent, which is the address in hex to int.
    addr = int(accessFields[2], 16)

    # derive the offset bits from address.
    memOffet = addr & (cacheLineSize - 1)

    # derive the index bits based on cache line size and # of sets.
    memSetIndex = addr >> int(math.log(cacheLineSize, 2)) & (numSets - 1)

    # derive tag by shifting bits.
    memTag = addr >> (int(math.log(numSets, 2)) + int(math.log(cacheLineSize, 2)))

    hitFlag = False

    set_start = linesPerSet * int(memSetIndex)
    set_end = set_start + linesPerSet

    for i in range(set_start, set_end):

        if (cache[i]["tag"] == memTag):
            hitFlag = True
            cache[i]["time"] = time.time()

            # write data to cache based on flag from trace file

            if(accessFields[1] == "W"):
                cache[i]["writeFlag"] = 1       # write back data from cache to main memory.

            break

    fullFlag = True

    min_time = 1000000000000000000000000  # use a large time value by default for LRU comparison
    LRU_slot = -1

    if(hitFlag == 0):
        missCount = missCount + 1
        for i in range(set_start, set_end):
            if (cache[i]["dataExists"] == 0):
                fullFlag = 0
                break
            if (cache[i]["time"] < min_time):
                min_time = cache[i]["time"]
                LRU_slot = i

        if (fullFlag == 0):
            cache[i]["tag"] = memTag
            cache[i]["dataExists"] = 1
            cache[i]["writeFlag"] = 0
            cache[i]["time"] = time.time()
        else:

            cache[LRU_slot]["tag"] = memTag
            cache[LRU_slot]["writeFlag"] = 0
            cache[LRU_slot]["dataExists"] = 1
            cache[LRU_slot]["time"] = time.time()

missRate = (missCount / hitCount) * 100
print("Cache miss rate: {:0.2f}%".format(missRate, 3))
