import numpy
from numpy.lib.type_check import isrealobj

def build_angle_dict(xtal_name, npz_dir):
    """
    Reads angle vectors I, J, K, theta for xtal_name from npz_dir
    Returns a dictionary allowing lookup of any existing bond angle ijk via angle_dict[i][j][k]
    """
    # the dictionary to return
    angle_dict = {}
    # read the 4 serialized vectors
    I = numpy.load(f"{npz_dir}/{xtal_name}_angles_I.npy")
    J = numpy.load(f"{npz_dir}/{xtal_name}_angles_J.npy")
    K = numpy.load(f"{npz_dir}/{xtal_name}_angles_K.npy")
    Theta = numpy.load(f"{npz_dir}/{xtal_name}_angles_theta.npy")
    # build the nested dictionaries
    for n in range(len(I)):
        if not I[n] in angle_dict.keys():
            angle_dict[I[n]] = {} # i level
        if not J[n] in angle_dict[I[n]].keys():
            angle_dict[I[n]][J[n]] = {} # j level
        angle_dict[I[n]][J[n]][K[n]] = Theta[n] # k level

    return angle_dict


def build_vector_dict(xtal_name, npz_dir):
    """
    Reads bond vector vectors I, J and matrix V for xtal_name from npz_dir
    Returns a dictionary allowing lookup of any existing bond vector ij via vector_dict[i][j]
    """
    vector_dict = {}

    I = numpy.load(f"{npz_dir}/{xtal_name}_vectors_I.npy")
    J = numpy.load(f"{npz_dir}/{xtal_name}_vectors_J.npy")
    V = numpy.load(f"{npz_dir}/{xtal_name}_vectors_V.npy")

    for n in range(len(I)):
        if not I[n] in vector_dict.keys():
            vector_dict[I[n]] = {}
        vector_dict[I[n]][J[n]] = V[:, n]

    return vector_dict


def build_distance_dict(xtal_name, npz_dir):
    """
    Reads vectors I, J, and D for xtal_name from npz_dir
    Returns a dictionary allowing lookup of any bond distance ij via distance_dict[i][j]
    """
    distance_dict = {}

    I = numpy.load(f"{npz_dir}/{xtal_name}_edges_src.npy")
    J = numpy.load(f"{npz_dir}/{xtal_name}_edges_dst.npy")
    D = numpy.load(f"{npz_dir}/{xtal_name}_euc.npy")

    for n in range(len(I)):
        if not I[n] in distance_dict.keys():
            distance_dict[I[n]] = {}
        distance_dict[I[n]][J[n]] = D[n]

    return distance_dict


if __name__ == "__main__":
    npz_dir = "data/graphs"
    xtal_name = "str_m7_o9_o24_bcu_sym.59"
    i = 148
    j = 62
    k = 78

    print(f"Loading data for {xtal_name} from {npz_dir}")

    distance_dict = build_distance_dict(xtal_name, npz_dir)
    angle_dict = build_angle_dict(xtal_name, npz_dir)
    vector_dict = build_vector_dict(xtal_name, npz_dir)

    #print(distance_dict)
    #print(angle_dict)
    #print(vector_dict)

    dij = distance_dict[i][j]
    dji = distance_dict[j][i]
    ijk = angle_dict[i][j][k]
    kji = angle_dict[k][j][i]
    aij = vector_dict[i][j]
    aji = vector_dict[j][i]

    print(f"i = {i}, j = {j}, k = {k}")
    print("Distance ij and distance ji: ", dij, dji)
    print("Angle ijk and angle kji: ",     ijk, kji)
    print("Vector ij and vector ji: ",     aij, aji)
