"""Program that takes in a time period defined in units used throughout the rest of the project and converts it to frequency units of 'per centimeter'"""

old_time=float(input("the time period you wish to convert is:"))

new_time=old_time*10.18e-15

freq=1/new_time

freq_centi=freq/2.99792458e10

print(freq_centi)
