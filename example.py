# import numpy as np
#
# dict = {'xMean': np.array([1])}
# print(type(dict['xMean']))
# print(np.concatenate((dict['xMean'], [2])) )
import numpy as np
a = np.array([1, 1, 1, 2, 2, 2, 2])
# print([item for item, count in collections.Counter(a).items() if count > 1])
# print([item for item in set(a) if a.count(item) > 1])
# b = list(set(a))
# print(b)
# unique, unique_indices, unique_inverse , unique_counts  = np.unique(a, return_index = True, return_inverse = True, return_counts = True)
# print(unique)
#
# print(unique_indices)
# print(unique_inverse)
# print(unique_counts)
#
# print(unique[unique_inverse])
# print(a[a == 1])
print(type([x**2 for x in a]))
