"""
Create serie of numbers and reverse it using FFT, conjugate and iFFT.
"""
import numpy as np
import matplotlib.pyplot as plt

# plot SX blue in complex plane and CSX red in separate function
def plot_complex_plane(SX, CSX):
    fig, ax = plt.subplots()
    ax.plot(SX.real, SX.imag, 'bo', alpha=0.5, label='SX')
    ax.plot(CSX.real, CSX.imag, 'ro', alpha=0.5, label='CSX')

    # plot a line between each point of SX and its corresponding point in CSX
    for sx, csx in zip(SX, CSX):
        ax.plot([sx.real, csx.real], [sx.imag, csx.imag], 'k--', alpha=0.5)

    ax.set_xlabel('Real')
    ax.set_ylabel('Imaginary')
    ax.legend()
    plt.show()

# set length of array
n = 16

# create random array with mean of 0
x = np.random.randint(-10, 10, size=n).astype('float64')

# subtract mean from all values
x -= np.mean(x) 

# round all values to integers
x = np.round(x).astype('int32')

# print x
print(x)

# calculate mean of x and print
print('mean: ', np.mean(x))

# calculate standard deviation of x
print('std: ', np.std(x))

# calculate variance of x and print
print('var: ', np.var(x))

# calculate fft of x
SX = np.fft.fft(x)

# compute complex conjugate of SX
CSX = np.conj(SX)

# print SX and CSX
#plot_complex_plane(SX, CSX)

# calculate ifft of CSX
iCSX = np.fft.ifft(CSX)

# reverse iCSX
riCSX = iCSX[::-1]

# print real part of iCSX and riCSX without decimals
print(np.real(iCSX).astype('int32'))
print(np.real(riCSX).astype('int32'))
