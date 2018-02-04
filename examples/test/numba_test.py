


if __name__ == '__main__':
    from numba import jit
    import numpy as np

    def sum2d(arr):
        M, N = arr.shape
        result = 0.0
        for i in range(M):
            for j in range(N):
                result += arr[i, j]
        return result

    print(sum2d(np.random.random_sample([10000,10000])))

