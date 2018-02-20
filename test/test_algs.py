import numpy as np
import smithWater  

def testSmithWaterman():
    """
    check if pos pair scores are greater than neg pair scores on average
    """
    assert smithWater.getAverageScores("MATIO", -7, -3)[0] > smithWater.getAverageScores("MATIO", -7, -3)[1]
    assert smithWater.getAverageScores("BLOSUM50", -7, -3)[0] > smithWater.getAverageScores("BLOSUM50", -7, -3)[1]

def testCombinations():
    """
    make sure I'm checking every combination for gap and extension penalty
    """
    assert len(smithWater.getCombinations()) == 100 

    
