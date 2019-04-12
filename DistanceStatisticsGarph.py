import os
import numpy as np
import defines
import GabiFuncs

import plotly.graph_objs as go
import plotly.offline as ply
import pickle



def plotGraph(distList):
    # Create sample data
    x = range(0, len(distList))
    y = distList
    # Create traces - data collections
    trace = go.Scatter(
        x = x,
        y = y,
        name = "Levenshtein Dist graph",
        line = dict (
            color = ("purple"),
            width = 4,
            # dash = 'dot'
        ),
    )
    # Create information / layout dictionary
    layout = dict (
        title = "Levenshtein Dist Graph of monitored data",
        xaxis = {'title' : 'bits shift in the monitored data' },
        yaxis = {'title' : 'Levenshtein Dist' }
    )
    # Pack the data
    data = [trace]
    # Create a figure
    fig = dict (data=data, layout=layout)
    # Plot
    ply.plot(fig, filename='./data/LevenshteinDistGraph.html')



# def hamming2(s1, s2):
#     """Calculate the Hamming distance between two bit strings"""
#     assert len(s1) == len(s2)
#     return sum(c1 != c2 for c1, c2 in zip(s1, s2))




def Compute_dist_by_shifts_bits_using_levenshtein(originalKeyStr, monitoredKeysStr):
    '''
    this function Compute dist by shifts bits using levenshtein in monitoredKeysStr
    
    :param originalKeyStr: the original key
    :param monitoredKeysStr: the monitored data in which we have to find string that matches 
        the original key with minimum levenshtein 
    :return listOfAllDicts: list of dicts with the distance, flips, deletions , insertion #{"DIST": _; "I": _; "F": _, "D"_X}
    :return distValuesList: list of all distance calculated
    :return numberOfAvalibleKeys: number of avalible keys  can be in the given monitoredKeysStr

    '''
    listOfAllDicts = []
    distValuesList = []
    
    #calculate how many keys can be in the given monitoredKeysStr
    # xxxxx =>5
    # yyyy =>4
    # 5-4=1
    #but there is 2 options so we add +1
    numberOfAvalibleKeys = len(monitoredKeysStr) - len(originalKeyStr) + 1
    originalKeyLen = len(originalKeyStr)
    for i in xrange(0, numberOfAvalibleKeys):
        print i
        tempStr = monitoredKeysStr[i: i+originalKeyLen]

        #dict = GabiFuncs.levenshtein_edit_dist(originalKeyStr, tempStr, show_strings=False)

        dict=({"DIST":1,"I":2,"D":3,"F":4},)
        if GabiFuncs.hamming_dist(originalKeyStr,tempStr) < len(originalKeyStr)/8:
            dict = GabiFuncs.levenshtein_edit_dist(originalKeyStr, tempStr, show_strings=False)
        else:
            dict[0]["DIST"]=len(originalKeyStr)/8
            dict[0]["I"] = len(originalKeyStr) / 8
            dict[0]["F"] = len(originalKeyStr) / 8
            dict[0]["D"] = len(originalKeyStr) / 8

        listOfAllDicts.append(dict[0])
        print dict[0]
        FILE_FOR_ALL_DICTS.write("%s\n" % dict[0])
        distValuesList.append(dict[0]["DIST"])
    return listOfAllDicts, distValuesList, numberOfAvalibleKeys




def FinedAllMinimumsDicts(listOfAllDicts, distValuesList, numberOfAvalibleKeys, len):
    """
        this function finds all the minimus distances in
         
    :param listOfAllDicts: 
    :param distValuesList: 
    :param numberOfAvalibleKeys: 
    :param len: 
    :return listOfAllDictsWithMinimumDist: 
    :return allMinmunsIndex: 
    """
    listOfAllDictsWithMinimumDist = []
    allMinmunsIndex = []
    distValuesListNumpy = np.array(distValuesList)
    firstMinIndex = distValuesListNumpy[0: 2*len].argmin()
    allMinmunsIndex.append(firstMinIndex)
    listOfAllDictsWithMinimumDist.append(listOfAllDicts[firstMinIndex])
    print listOfAllDicts[firstMinIndex]
    firstMinIndex = firstMinIndex + len/2
    stratIndex = firstMinIndex
    while stratIndex <= numberOfAvalibleKeys:
        endIndex = stratIndex + len + 1
        if endIndex <= numberOfAvalibleKeys:
            reletiveIndex = distValuesListNumpy[stratIndex: endIndex].argmin()
            globalIndex = stratIndex + reletiveIndex
            allMinmunsIndex.append(globalIndex)
            listOfAllDictsWithMinimumDist.append(listOfAllDicts[globalIndex])
            FILE_FOR_ALL_MINIMUM_DICTS.write("%s\n" % listOfAllDicts[globalIndex])
            print listOfAllDicts[globalIndex]
            stratIndex = endIndex
        else:
            return listOfAllDictsWithMinimumDist, allMinmunsIndex
    return listOfAllDictsWithMinimumDist, allMinmunsIndex





def CaclulateAvg(listOfAllDictsWithMinimumDist):
    '''
    
    :param listOfAllDictsWithMinimumDist: list of dictionaries that holds computeation of all the minimum errors in monitored data string collected from cache LLC side channel using levenshtein distance
        (error rate includes distance, flips, deletions , insertion) 
    :return: dict with the  average distance, average flips, average deletions , average insertion
        {"insertionAVG": _; "deletionAVG": _; "flipsAVG": _, "LevDistAVG" _}
    '''
    insertionList = list(map(lambda x: x['I'], listOfAllDictsWithMinimumDist))
    deletionList = list(map(lambda x: x['D'], listOfAllDictsWithMinimumDist))
    flipsList = list(map(lambda x: x['F'], listOfAllDictsWithMinimumDist))
    distsList = list(map(lambda x: x['F'], listOfAllDictsWithMinimumDist))

    insertionListNumpy = np.array(insertionList)
    deletionListNumpy = np.array(deletionList)
    flipsListNumpy = np.array(flipsList)
    distsListNumpy = np.array(distsList)

    result = {}
    result["insertionAVG"] = np.mean(insertionListNumpy)
    result["deletionAVG"] = np.mean(deletionListNumpy)
    result["flipsAVG"] = np.mean(flipsListNumpy)
    result["LevDistAVG"] = np.mean(distsListNumpy)
    return result



def ComputeAvgOfAllMinumumsInMonitoredKeys(originalKeyStr, monitoredKeysStr):
    '''
    This is the main function which computes 
    the average error in monitored data string collected from cache LLC side channel using levenshtein distance
    error rate includes average distance, average flips, average deletions , average insertion 
    and plots graph 
    :param originalKeyStr: the original key
    :param monitoredKeysStr: the monitored data in which we have to find string that matches 
        the original key with minimum levenshtein 
    :return: None
    '''
    
    conclusiton = {}

    listOfAllDicts, distValuesList, numberOfAvalibleKeys = Compute_dist_by_shifts_bits_using_levenshtein(originalKeyStr, monitoredKeysStr)

    str = "numberOfAvalibleKeys ="
    print str
    FILE_FOR_MORE_INFORMATION.write("%s\n" % str)
    print numberOfAvalibleKeys
    FILE_FOR_MORE_INFORMATION.write("%s\n" % numberOfAvalibleKeys)

    plotGraph(distValuesList)
    listOfAllDictsWithMinimumDist, allMinmunsIndex = FinedAllMinimumsDicts(listOfAllDicts, distValuesList, numberOfAvalibleKeys, len(originalKeyStr))

    str = "the indexs of all minimus values are:"
    print str
    FILE_FOR_MORE_INFORMATION.write("%s\n" % str)
    print allMinmunsIndex
    FILE_FOR_MORE_INFORMATION.write("%s\n" % allMinmunsIndex)

    conclusiton = CaclulateAvg(listOfAllDictsWithMinimumDist) 

    str = "THE AVRAGE OF NOISE IS:"
    print str
    FILE_FOR_MORE_INFORMATION.write("%s\n" % str)
    print conclusiton
    FILE_FOR_MORE_INFORMATION.write("%s\n" % conclusiton)
    return conclusiton




if __name__ == "main":
    FILE_FOR_ALL_DICTS = open('./data/listOfAllDicts.txt', 'w')
    FILE_FOR_ALL_MINIMUM_DICTS = open('./data/listOfAllDictsWithMinimumDist.txt', 'w')
    FILE_FOR_MORE_INFORMATION = open('./data/moreInfo.txt', 'w')
        
    originalKeyFile = open(defines.FILE_ORIGINAL_KEY, "rt") #open file
    # originalKeyStr = originalKeyFile.readline()
    originalKeyStr = originalKeyFile.read().replace('\n','').replace(' ','').replace(',','')
    originalKeyFile.close()
    str = "originalKeyStr len= "
    print str
    FILE_FOR_MORE_INFORMATION.write("%s\n" % str)
    print len(originalKeyStr)
    FILE_FOR_MORE_INFORMATION.write("%s\n" % len(originalKeyStr))
    
    
    monitoredaKeysFile = open(defines.FILE_SAMPLED_DATA, "rt") #open file
    # monitoredKeysStr = monitoredaKeysFile.readline()
    monitoredKeysStr = monitoredaKeysFile.read().replace('\n','').replace(' ','').replace('[','').replace(']','').replace(',','')
    monitoredaKeysFile.close()
    
    str = "monitoredKeysStr len= "
    print str
    FILE_FOR_MORE_INFORMATION.write("%s\n" % str)
    print len(monitoredKeysStr)
    FILE_FOR_MORE_INFORMATION.write("%d\n" % (len(monitoredKeysStr)))
    
    
    monitoredKeysStr = monitoredKeysStr[0:len(monitoredKeysStr)/1]
    str = "new monitoredKeysStr len= "
    print str
    FILE_FOR_MORE_INFORMATION.write("%s\n" % str)
    print len(monitoredKeysStr)
    FILE_FOR_MORE_INFORMATION.write("%d\n" % (len(monitoredKeysStr)/1))
    
    
    result = ComputeAvgOfAllMinumumsInMonitoredKeys(originalKeyStr, monitoredKeysStr)
    
    FILE_FOR_ALL_DICTS.close()
    FILE_FOR_ALL_MINIMUM_DICTS.close()
    FILE_FOR_MORE_INFORMATION.close()





