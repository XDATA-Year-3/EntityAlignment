import scipy.io as sio
# Reads a sparse Matlab matrix stored in a .mat file and returns the non-zero entries as a
# dictionary with vertex number from graph A as the key and a list of matching vertices
# along with weight as the value
#
# INPUTS :  matFilename : path to the .mat file to be read
#               varName : name of the matrix that was stored in the .mat file
# OUTPUTS :     matches : a dictionary with the vertex index from graph A as the key
#                         and the value is a list of tuples where the second entry of the
#                         tuple is the potential match from graph B and the first entry
#                         is the weight of this match
def readTopK(matFilename, varName):
    wkspace = sio.loadmat(matFilename);
    # Because Matlab stores the entire workspace in the .mat file, it is necessary to
    # know the name used in Matlab to store the matrix.
    # Is there a better way to do this or a workaround?
    topkMatrix = wkspace[varName];

    matches = {}; # dictionary of matches

    # lists of corresponding indices of nonzero entries (i.e., potential matches!)
    nonZeroA, nonZeroB = topkMatrix.nonzero();

    for i,j in zip(nonZeroA, nonZeroB): # zip is really useful
        if (i in matches):
            matches[i].append((topkMatrix[i,j], j))
        else:
            matches[i] = [(topkMatrix[i,j], j)];

    return matches

