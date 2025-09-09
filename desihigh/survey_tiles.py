import numpy as np

def generate_tile_data():

    # right ascension coordinates of tiles for the 2021-09-22 DESI observing plan
    ra = [
       266.0, 261.4, 312.4, 317.0, 318.7, 321.7, 326.3, 327.6, 329.6, 333.0, 336.3, 337.7, 340.1, 
       337.4, 338.6, 336.3, 335.8, 336.1, 335.2, 335.3, 336.6, 96.6, 94.3, 98.5, 102.0,
       103.2, 103.8, 105.3, 106.6, 109.1, 109.3
    ]
    
    # transform the right ascention values to fall between 110 and -100 degrees 
    ra = np.array(ra)
    ra = (ra - 150)%360 + 150 - 360
    
    # declination coordinates of tiles for the 2021-09-22 DESI observing plan
    declination = [
       24.8, 12.8, 0.5, -2.7, 2.3, 4.1, 0.6, -2.5, 0.2, -0.2, -0.6, -6.4, -12.6, 26.1, 31.9, 19.4, 
        15.3, 29.3, 23.4, 9.3, 6.2, 62.6, 65.9, 69.0, 64.9, 61.3, 52.7, 49.5, 43.7 , 35.3, 40.0
    ]
    
    ra.tofile('../data/20210922_tiles_ra.BIN')
    np.array(declination).tofile('../data/20210922_tiles_dec.BIN')
    
