#!/opt/python2.7/bin/python
# above on bhmgiap09lt
## on bhmgiapp14ld: /usr/bin/python

import numpy as np
import pandas as pd

# from https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe/41078786#41078786
# Sometimes you have text based data in there as well. This lets you specify 
# the columns cols of your numerical data and will run quantile normalization 
# on those columns. At the end it will merge back the non-numeric (or not to be
# normalized) columns from your original data frame. 
#
# copy dataframe and only use the columns with numerical values
# copy() copies data from inputs. Only affects DataFrame / 2d ndarray input
#
# filter():
# Subset rows or columns of dataframe according to labels in the
#    specified index.
# Note that this routine does not filter a dataframe on its contents.
#    The filter is applied to the labels of the index.
# items=cols: List of info axis to restrict to (must not all be present)
#   List of info axis to restrict to (must not all be present). Axis default
#   is columns

def qn(dataframe, cols):
    print 'quantileNormalize.py dataframe.to_dict(): %s' % dataframe.to_dict()
    print 'quantileNormalize.py cols: %s' % cols

    df = dataframe.copy().filter(items=cols)

    # columns from the original dataframe not specified in cols
    non_numeric = dataframe.filter(items=list(filter(lambda col: col not in cols, list(dataframe))))

    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.stack.html
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.groupby.html
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html
     # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.rank.html

    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()  

    norm = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()


    result = pd.concat([norm, non_numeric], axis=1)
    return result
