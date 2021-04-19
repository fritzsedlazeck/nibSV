# vim: sts=4:ts=4:sw=4:et:tw=0
#from cpuinfo import nil
from math import nil
from os import nil
#from threadpool import nil
from streams import nil
from strformat import fmt
from strutils import nil
import heapqueue
import osproc
import times

type PbError* = object of CatchableError
type GenomeCoverageError* = object of PbError
type FieldTooLongError* = object of PbError
type TooFewFieldsError* = object of PbError

proc raiseEx*(msg: string) {.discardable.} =
    raise newException(PbError, msg)

proc isEmptyFile*(fn: string): bool =
    var finfo = os.getFileInfo(fn)
    if finfo.size == 0:
        return true
    return false

#from strformat import fmt
proc isOlderFile*(afn, bfn: string): bool =
    ## Return true iff afn is older than bnf.
    let
        at = os.getLastModificationTime(afn)
        bt = os.getLastModificationTime(bfn)
        #af = at.format("yyyy-MM-dd'T'HH:mm:ss,ffffffzzz")
        #bf = bt.format("yyyy-MM-dd'T'HH:mm:ss,ffffffzzz")
    #echo "glmt {afn}: {af}, {bfn}: {bf}".fmt
    return at < bt

template withcd*(newdir: string, statements: untyped) =
    let olddir = os.getCurrentDir()
    os.setCurrentDir(newdir)
    defer: os.setCurrentDir(olddir)
    statements

proc log*(words: varargs[string, `$`]) =
    for word in words:
        write(stderr, word)
    write(stderr, '\l')

proc logt*(words: varargs[string, `$`]) =
    var then {.global.} = times.now()
    let
        since = times.initDuration(seconds = times.inSeconds(times.now() - then))
        dp = times.toParts(since)
        prefix = strformat.fmt("{dp[Hours]}:{dp[Minutes]:02d}:{dp[Seconds]:02d}s ")
    write(stderr, prefix)
    log(words)

proc adjustThreadPool*(n: int) =
    ## n==0 => use ncpus
    ## n==-1 => do not alter threadpool size (to avoid a weird problem for now)
    log("(ThreadPool is currently not used.)")
    #var size = n
    #if n == 0:
    #    size = cpuinfo.countProcessors()
    #if size > threadpool.MaxThreadPoolSize:
    #    size = threadpool.MaxThreadPoolSize
    #if size == -1:
    #    log("ThreadPoolsize=", size,
    #        " (i.e. do not change)",
    #        ", MaxThreadPoolSize=", threadpool.MaxThreadPoolSize,
    #        ", NumCpus=", cpuinfo.countProcessors())
    #    return
    #log("ThreadPoolsize=", size,
    #    ", MaxThreadPoolSize=", threadpool.MaxThreadPoolSize,
    #    ", NumCpus=", cpuinfo.countProcessors())
    #threadpool.setMaxPoolSize(size)

iterator walk*(dir: string, followlinks = false, relative = false): string =
    ## similar to python os.walk(), but always topdown and no "onerror"
    # Slow! 30x slower than Unix find.
    let followFilter = if followLinks: {os.pcDir, os.pcLinkToDir} else: {os.pcDir}
    let yieldFilter = {os.pcFile, os.pcLinkToFile}
    for p in os.walkDirRec(dir, yieldFilter = yieldFilter,
            followFilter = followFilter, relative = relative):
        yield p

iterator readProc*(cmd: string): string =
    ## Stream from Unix subprocess, e.g. "find .".
    ## But if cmd=="-", stream directly from stdin.
    if cmd == "-":
        log("Reading from stdin...")
        for line in lines(stdin):
            yield line
    else:
        log("Reading from '" & cmd & "'...")
        var p = osproc.startProcess(cmd, options = {poEvalCommand})
        if osproc.peekExitCode(p) > 0:
            let msg = "Immedate failure in readProc startProcess('" & cmd & "')"
            raiseEx(msg)
        defer: osproc.close(p)
        for line in streams.lines(osproc.outputStream(p)):
            yield line

iterator readProcInMemory(cmd: string): string =
    ## Read from Unix subprocess, e.g. "find .", into memory.
    ## But if cmd=="-", stream directly from stdin.
    if cmd == "-":
        log("Reading from stdin...")
        for line in lines(stdin):
            yield line
    else:
        log("Reading from '" & cmd & "'...")
        let found = osproc.execProcess(cmd, options = {poEvalCommand})
        var sin = streams.newStringStream(found)
        for line in streams.lines(sin):
            yield line

proc removeFile*(fn: string, failIfMissing = false) =
    if failIfMissing and not os.fileExists(fn):
        raiseEx("Cannot remove non-existent file '" & fn & "'")
    log("rm -f ", fn)
    os.removeFile(fn)

proc removeFiles*(fns: openarray[string], failIfMissing = false) =
    for fn in fns:
        removeFile(fn, failIfMissing)

proc which*(exe: string) =
    let cmd = "which " & exe
    log(cmd)
    discard execCmd(cmd)

proc thousands*(v: SomeInteger): string =
    if v == 0:
        return "0"
    var i: type(v) = v
    let negative = (i < 0)
    i = abs(i)
    #result = strformat.fmt"{i mod 1000:03}"
    #i = i div 1000
    while i > 0:
        result = strformat.fmt"{i mod 1000:03}," & result
        i = i div 1000
    # Drop tailing comma.
    assert result[^1] == ','
    result = result[0 .. ^2]
    # Drop leading 0s.
    while result[0] == '0':
        result = result[1 .. ^1]
    if negative:
        result = '-' & result

proc splitWeighted*(n: int, sizes: seq[int]): seq[int] =
    # Split sizes into n contiguous subsets, weighted by each size.
    # Each elem of result will represent a range of elems of sizes.
    # len(result) will be <= n

    if n == 0:
        return
    var sums: seq[int]
    var totalSize = math.sum(sizes)
    var remSize = totalSize
    var curr = 0
    var remN = min(n, len(sizes))
    while len(sizes) > curr:
        #assert len(sizes) > curr, "not enough elements in sizes {len(sizes)} <= {curr}".fmt
        result.add(0)
        let approx = int(math.ceil(remSize / remN))
        #echo "approx={approx}, remaining={remN}, tot={remSize}".fmt
        sums.add(0)
        while sums[^1] < approx:
            result[^1] += 1
            sums[^1] += sizes[curr]
            curr += 1
        remN -= 1
        remSize -= sums[^1]
    assert math.sum(result) == len(sizes)
    assert math.sum(sizes) == totalSize
    assert len(result) <= n

type
    BinSum = object
        indices: seq[int]
        sum: int64
        order: int
    WeightedIndex = tuple[index: int, size: int]

proc `<`(a, b: BinSum): bool =
    return a.sum < b.sum or (a.sum == b.sum and a.indices.len() < b.indices.len()) or
        (a.sum == b.sum and a.indices.len() == b.indices.len() and a.order < b.order)
proc `<`(a, b: WeightedIndex): bool =
    return a.size > b.size or (a.size == b.size and a.index > b.index)

proc partitionWeighted*(n: int, sizes: seq[int]): seq[seq[int]] =
    ## {sizes} is an index; other seqs refer to its indices.
    ## The splits for this version are not required to be contiguous.
    ## The result has at most n index-seqs, none of which are empty.
    var biggest = initHeapQueue[WeightedIndex]()
    for i in 0 ..< len(sizes):
        let wi: WeightedIndex = (index: i, size: sizes[i])
        biggest.push(wi)
    var smallest_bin = initHeapQueue[BinSum]()
    for x in 0 ..< n:
        var bin: BinSum = BinSum(sum: 0, order: x)
        smallest_bin.push(bin)
    while biggest.len() > 0:
        let wi = biggest.pop()
        var bin = smallest_bin.pop()
        bin.indices.add(wi.index)
        bin.sum += wi.size
        smallest_bin.push(bin)
    while smallest_bin.len() > 0:
        let bin = smallest_bin.pop()
        if bin.indices.len() > 0:
            result.add(bin.indices)
    return result

proc combineToTarget*(target: int64, weights: seq[int64]): seq[seq[int]] =
    # Given a seq of weights,
    # combine consecutive groups of them until they meet target.
    # Return a seq of seqs of those indices. For now,
    # the results will always be consecutive, e.g.
    # [ [0,1,2], [2,3], [4], [5,6] ]
    var
        total = target
        n = -1
    for i in 0 ..< len(weights):
        let next_weight = weights[i]
        #echo "i:{i} next:{next_weight} total:{total} n:{n}".fmt
        if total >= target:
            # new group
            result.add(@[i])
            n = len(result) - 1
            total = next_weight
        else:
            # current group
            result[n].add(i)
            total += next_weight

const
    MAX_HEADROOM* = 1024
type
    Headroom* = array[MAX_HEADROOM, cchar]

proc sscanf*(s: cstring, frmt: cstring): cint {.varargs, importc,
        header: "<stdio.h>".}

proc strlen(s: cstring): cint {.importc: "strlen", nodecl.}

proc strlen(a: var Headroom): int =
    let n = strlen(cast[cstring](addr a))
    return n

proc toString*(ins: var Headroom, outs: var string, source: string = "") =
    var n = strlen(ins)
    if n >= (MAX_HEADROOM - 1):
        # Why is max-1 illegal? B/c this is used after sscanf, and that has no way to report
        # a buffer-overflow. So a 0 at end-of-buffer is considered too long.
        let msg = strformat.fmt"Too many characters in substring (>{MAX_HEADROOM - 1}) from '{source}'"
        raise newException(util.FieldTooLongError, msg)
    outs.setLen(n)
    for i in 0 ..< n:
        outs[i] = ins[i]

proc getNthWord*(line: string, n: Natural, delim: char): string =
    ## n is 0-based
    var
        start = 0
        count = 0
        found = -1
    while count < n:
        found = strutils.find(line, delim, start)
        if found == -1:
            let msg = "Found only {count} < {n} instances of '{delim}' in '{line}'".fmt
            raiseEx(msg)
        start = found + 1
        count += 1
    var wordEnd = strutils.find(line, delim, start)
    if wordEnd == -1:
        wordEnd = line.len()
    return line[start..(wordEnd-1)]
