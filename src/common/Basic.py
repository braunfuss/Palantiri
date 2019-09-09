import os
import sys
import time
import shutil
import subprocess
import getpass
from palantiri.common import Logfile
import httplib2

if sys.version_info.major >= 3:
    from urllib.request import urlopen
    import _pickle as pickle
else:
    from urllib2 import urlopen
    import cPickle as pickle


def floatToString(fList, format=None, delim=','):
    if not format:
        return delim.join(map(str, fList))
    else:
        s = []

        for val in fList:
            s.append(format % val)
        return delim.join(s)


def stringToFloat(s, delim=','):

    words = s.split(delim)
    line = []

    for i in range(len(words)):
        if words[i] == '\n':
            break

        line.append(float(words[i]))

    return line


def matrixToString(matrix, nLines, nColumns, format=None, delim=','):

    lines = []

    for i in range(nLines):
        lines.append(floatToString(matrix[i], format))

    return '\n'.join(lines)


def stringToMatrix(lines, nLines, nColumns, delim=','):

    matrix = []

    for i in range(nLines):
        vect = stringToFloat(lines[i], delim)
        assert len(vect) == nColumns
        matrix.append(vect)

    assert len(matrix) == nLines
    return matrix


def  writeVector(fileName, vector, format=None):
        writeTextFile(fileName, list(floatToString(vector, format)))


def readVector(fileName):
    return stringToFloat(readTextFile(fileName, 1)[0])


def writeMatrix(fileName, matrix, nLines, nColumns, format=None):
    writeTextFile(fileName, matrixToString(matrix, nLines, nColumns, format))


def readMatrix(fileName, nLines, nColumns):
    lines = readTextFile(fileName)
    matrix = stringToMatrix(lines,  nLines, nColumns)

    return matrix


def formatStrings(strings, format1):

    result = []

    try:
        for i in range(len(strings)):
            result.append(format1 % (strings[i]))

    except Exception:
        Logfile.exception('formatStrings', 'Illegal format', abortProg=True)

    return result


def selectStrings(strings, mask):

    result = []

    for i in range(len(mask)):
        if mask[i]:
            result.append(strings[i])

    return result


def _stringsEndsWith(strings, postfixList):

    assert len(postfixList) > 0
    mask = []

    for i in range(len(strings)):
        s = strings[i]
        b = False

        for postfix in postfixList:
            if s.endswith(postfix):
                b = True
                break

        mask.append(b)

    assert len(mask) == len(strings)
    return mask


def stringsEndsWith(strings, postfixList):

    if type(postfixList) is str:
        list1 = list(postfixList)
    else:
        list1 = postfixList

    return _stringsEndsWith(strings, list1)


def toStringList(arg0, arg1=None, arg2=None, arg3=None, arg4=None):

    s = []
    s.append(arg0)

    if arg1 is not None:
        s.append(arg1)
    if arg2 is not None:
        s.append(arg2)
    if arg3 is not None:
        s.append(arg3)
    if arg4 is not None:
        s.append(arg4)

    return s


def Not(mask):

    result = []

    for i in range(len(mask)):
        if mask[i]:
            result.append(False)
        else:
            result.append(True)

    return result


def And(mask):

    for i in range(len(mask)):
        if not mask[i]:
            return False

    return True


def baseFileName(fullName):
    basename = os.path.basename(fullName)
    filename = os.path.splitext(basename)
    return filename[0]


def isNumber(s):

    try:
        float(s)
    except:
        return False

    return True


def isInt(s):

    if not isNumber(s):
        return False
    try:
        int(s)
    except:
        return False

    return True


def checkIsNumber(string, minVal=None, maxVal=None):

    assert isNumber(minVal) and isNumber(maxVal)

    msg = None
    s1 = 'Value ' + string + ' '

    if not isNumber(string):
        msg = s1 + 'is not a number'
    else:
        val = float(string)

        if minVal is None and maxVal is None:
            msg = None
        elif minVal is None and val > maxVal:
            msg = s1 + '> ' + str(maxVal)
        elif maxVal is None and val < minVal:
            msg = s1 + '< ' + str(minVal)

        elif val < minVal or val > maxVal:
            msg = s1 + 'outside range [' + str(minVal) + ',' +\
                     str(maxVal) + ']'
    return msg


def checkGreaterZero(string):

    msg = None
    s1 = 'Value ' + string + ' '

    if not isNumber(string):
        msg = s1 + 'is not a number'
    else:
        val = float(string)

        if val is 0.0:
            msg = s1 + 'is zero'
        elif val < 0.0:
            msg = s1 + '< 0.0'
    return msg


def checkNotNegative(string):

    msg = None
    s1 = 'Value ' + string + ' '

    if not isNumber(string):
        msg = s1 + 'is not a number'
    elif float(string) < 0.0:
        msg = s1 + '< 0.0'

    return msg


def checkExistsKeys(dict, keyList, isAbort=False):

    isOk = True

    for key in keyList:
        if not key in dict:
            isOk = Logfile.error('Key <' + str(key) + '> missing in config file')

    if isOk:
        return True
    elif isAbort:
        Logfile.abort()

    return False


def checkExistsDir(dirName, isAbort=False):

    if os.path.isdir(dirName):
        return True

    Logfile.error('Cannot find directory', dirName)

    if isAbort:
        Logfile.abort()
    return False


def createDirectory(dirName, optional=False):

    if os.path.isdir(dirName):
        return True

    os.makedirs(dirName)

    if os.path.isdir(dirName):
        return True
    else:
        Logfile.error('Cannot open directory', dirName)

    if not optional:
        Logfile.abort()

    return False


def changeDirectory(dirPath):

    path = []

    if not type(dirPath) is list:
        path.append(dirPath)
    else:
        path = dirPath

    for dir in path:
        createDirectory(dir)
        os.chdir(dir)

    return os.getcwd()


def checkFileExists(fileName, isAbort=False):

    if os.path.isfile(fileName):
        return True

    Logfile.fileOpenError(fileName)

    if isAbort:
        Logfile.abort('')
    else:
        return False


def openTextFile(fileName, mode):

    if mode == 'r' and not checkFileExists(fileName):
        return None

    return open(fileName, mode)


def readTextFile(fileName, maxLines=-1):

    lines = []
    fp = openTextFile(fileName, 'r')

    if fp is None:
        return lines

    if maxLines is -1:
        lines = fp.readlines()
    else:
        lines = fp.readlines(maxLines)

    fp.close()

    return lines


def writeTextFile(fileName, lines):

    fp = openTextFile(fileName, 'w')
    if fp is None:
        return

    for s in lines:
        fp.write(s)

    fp.close()


def appendToFile(fileName, lines):

    fp = openTextFile(fileName, 'a')
    if fp is None:
        return

    for s in lines:
        fp.write(s)

    fp.close()


def copyFile(srcFile, destFile, isAbort=False):

    if not checkFileExists(srcFile, isAbort):
        return False

    shutil.copy(srcFile, destFile)
    return True


def dumpToFile(fileName, anyData):
    if sys.version_info.major >= 3:
        output = open(fileName, 'wb')
    else:
        output = open(fileName, 'wb')

    pickle.dump(anyData, output)
    output.close()


def loadDump(fileName):

    if not os.path.isfile(fileName):
        Logfile.fileOpenError(fileName)
        return None

    if sys.version_info.major >= 3:
        f = open(fileName, 'rb')
    else:
        f = open(fileName)
    data = pickle.load(f)
    f.close()
    os.remove(fileName)

    return data


def removeFile(files):
    if not files is list:
        files2 = [files]
    else:
        files2 = files

    for f in files2:
        if os.path.isfile(f):
            os.remove(f)


def removeFiles(dir, prefix=None):

    names = os.listdir(dir)
    names2 = []

    if prefix is None:
        names2 = names
    else:
        for file in names:
            if file.startswith(prefix):
                names2.append(file)

    for file in names2:
        if dir != '.':
            fullName = os.path.join(dir, file)
        else:
            fullName = file

        if os.path.isfile(fullName):
            if dir != '.':
                file1 = dir + '/' + file
            else:
                file1 = file

            os.remove(file1)


def wait(seconds, prompt=True):

    for i in range(seconds):
        if prompt:
            print(str(i))
        time.sleep(1)


def systemCmd(cmd):

    if type(cmd) is list:
        cmd2 = ' '.join(cmd)
    else:
        cmd2 = cmd

    pipe = subprocess.Popen(cmd2, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT).stdout
    lines = pipe.readlines()
    lines2 = []

    for s in lines:
        if s.endswith('\r\n'):
            lines2.append(s[:-2])
        else:
            lines2.append(s)

    return lines2


def killByPID(pidList):

    ownPID = os.getpid()

    for pid in pidList:
        if pid == ownPID:
            continue

        os.system('kill ' + str(pid))


def removeTempFiles():

    dir = '/tmp'
    names = os.listdir(dir)
    user = getpass.getuser()
    cnt = 0

    for s in names:
        if s .startswith('obspy-'):
            os.remove(dir + '/' + s)
            cnt += 1

    if cnt != 0:
        Logfile.add('Remove ' + str(cnt) + ' temp files obspy-*.*')

    return cnt


def existsHTML_Page(url, text=None, withComment=False, isAbort=False):


    h = httplib2.Http()

    try:
        resp = h.request(url, 'HEAD')
        status = int(resp[0]['status'])

        if status < 400:
            return True   # page exists
    except:
        status = 104  # error: [Errno 104] Connection reset by peer  ???

    if withComment:
        s = 'HTML page not found, error = ' + str(status)

        if text is None:
            Logfile.error(s, url)
        else:
            Logfile.error(text, s, url)

    if isAbort:
        Logfile.abort()

    return False


def readURL(url, tmpFile=None):

    lines = []

    try:
        datasource = urlopen(url)

        while True:
            line = datasource.readline()

            if line == '':
                break
            lines.append(line)

        if len(lines) == 0:
            Logfile.error('Cannot read from URL:', url)
            lines = []

        Logfile.add(str(len(lines)) + ' lines read from URL:', url)

        if tmpFile:
            writeTextFile(tmpFile, lines)

    except:
        Logfile.exception('readURL')

    return lines


def readUrl2(pythonScript, pythonScript_2):

    assert pythonScript.endswith('.py')

    try:
        lines = []
        url = 'http://www.staedtke-berlin.kilu.de/' + pythonScript
        page = urlopen(url)

        for line in page:
            lines.append(line[:-1])

    except:
        return False

    try:
        writeTextFile(pythonScript_2, lines)
    except:
        return False

    return True


def question(text):

    while True:
        if sys.version_info.major >= 3:
            c = input(text + ' <y/n>:')
            c = c.lower()

            print('c=<' + c + '>')

            if c == 'y':
                return True
            elif c == 'n':
                return False
            else:
                continue
        else:
            c = raw_input(text + ' <y/n>:')
            c = c.lower()

            print('c=<' + c + '>')

            if c == 'y':
                return True
            elif c == 'n':
                return False
            else:
                continue


def getUserName():

    s = os.environ['USERDOMAIN']

    return s
