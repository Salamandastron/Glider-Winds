import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

#Create arrays to store the data in
time_array = []
lat_array = []
long_array = []
alt_pres_array = []
alt_gps_array = []

def rad_deg(theta):
    phi = theta*(180/(np.pi))
    return phi

def deg_rad(phi):
    theta = phi*((np.pi)/180)
    return theta

def difference(array):
    delta = array[i+1] - array[i]
    return delta


with open("JohnEllis_17-09-19.igc", "r") as fin:
    with open("John.txt", "w") as fout:
        for line in fin:
            if line.startswith('B'):
                time_hour = int(line[1:3])  #Extract the hour
                time_min = int(line[3:5])   #Extract the minute
                time_sec = int(line[5:7])   #Extract the second

                time = (time_hour*60*60) + (time_min*60) + time_sec #Seconds since midnight
                time_array.append(time) #Put seconds since midnight into an array
               
                if line[14:15] == 'S': #North is +ve south is -ve
                    lat_deg = int(line[7:9])*-1 + (int(line[9:14])*-1)/60000
                    "lat_min = (int(str[9:14])*-1)/60000"
                else:
                    lat_deg = int(line[7:9]) + int(line[9:14])/60000
                    "lat_min = int(str[9:14])/60000"

                lat_rad = deg_rad(lat_deg) #Convert to radians
                lat_array.append(lat_rad)
               
                if line[23:24] == 'W': #East is +ve west is -ve
                    long_deg = int(line[15:18])*-1 + (int(line[18:23])*-1)/60000
                    "long_min = (int(str[18:23])*-1)/60000"
                else:
                    long_deg = int(line[15:18]) + int(line[18:23])/60000
                    "long_min = int(str[18:23])/60000"

                long_rad = deg_rad(long_deg) #Convert to radians
                long_array.append(long_rad)
               
                alt_pres = int(line[25:30])
                alt_pres_array.append(alt_pres)
               
                alt_gps = int(line[30:35])
                alt_gps_array.append(alt_gps)
               
                fout.write(line)
            else:
                fout.write('\n')
    fout.close()
fin.close()

#Return the coordinates of a point clicked on by the the user on a graph
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    rx = round(ix)
    print ('x =', rx, 'y =', iy)

    global coords
    coords = [rx, iy]

    plt.close()

    return coords

#Binary search to find the index of a value within an array
#Also returns the index of the value closest (rounded up) to the inupt value if the input value is not in the array
def find_index(array, value):
    left, right = 0, len(array) -1

    while left <= right:
        middle = (left + right)//2

        if array[middle] == value:
            return middle
    
        if array[middle] < value:
            left = middle + 1

        if array[middle] > value:
            right = middle - 1

    return middle

#Circle fitting itterative method, known as Taubin-Newton  
def Taubin_Newton (x_pos, y_pos):
    nPoints = len(x_pos)
    if (nPoints < 3):
        print("Too few points")
    else:
        x_centroid = np.mean(x_pos)
        y_centroid = np.mean(y_pos)
        Mxx = 0 
        Myy = 0
        Mxy = 0
        Mxz = 0
        Myz = 0
        Mzz = 0
        i = 0
        while (i < nPoints):
            Xi = x_pos[i] - x_centroid
            Yi = y_pos[i] - y_centroid
            Zi = Xi * Xi + Yi * Yi
            Mxy += Xi * Yi
            Mxx += Xi * Xi
            Myy += Yi * Yi
            Mxz += Xi * Zi
            Myz += Yi * Zi
            Mzz += Zi * Zi
            i = i+1

        Mxx /= nPoints
        Myy /= nPoints
        Mxy /= nPoints
        Mxz /= nPoints
        Myz /= nPoints
        Mzz /= nPoints

        Mz = Mxx + Myy
        Cov_xy = Mxx * Myy - Mxy * Mxy
        A3 = 4 * Mz
        A2 = -3 * Mz * Mz - Mzz
        A1 = Mzz * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz - Mz * Mz * Mz
        A0 = Mxz * Mxz * Myy + Myz * Myz * Mxx - Mzz * Cov_xy - 2 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy
        A22 = A2 + A2
        A33 = A3 + A3 + A3

        xnew = 0
        ynew = 1e+20
        epsilon = 1e-12
        iterMax = 20

        iter = 0
        while (iter < iterMax):
            yold = ynew
            ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3))
            if (np.abs(ynew) > np.abs(yold)):
                System.out.println("Newton-Taubin goes wrong direction: |ynew| > |yold|")
                xnew = 0
                break
            
            Dy = A1 + xnew * (A22 + xnew * A33)
            xold = xnew
            xnew = xold - ynew / Dy
            if (np.abs((xnew - xold) / xnew) < epsilon):
                break
            
            if (iter >= iterMax):
                System.out.println("Newton-Taubin will not converge")
                xnew = 0
            
            if (xnew < 0.):
                System.out.println("Newton-Taubin negative root: x = " + xnew)
                xnew = 0

            iter = iter + 1
                
        centreRadius = [0, 0, 0, 0]
        det = xnew * xnew - xnew * Mz + Cov_xy
        x = (Mxz * (Myy - xnew) - Myz * Mxy) / (det * 2)
        y = (Myz * (Mxx - xnew) - Mxz * Mxy) / (det * 2)
        centreRadius[0] = x + x_centroid
        centreRadius[1] = y + y_centroid
        centreRadius[3] = np.sqrt(centreRadius[0]**2 + centreRadius[1]**2)

        return centreRadius

#Create and array for the bearings
bearing = [0] 
i = 0
#Calculate the bearings
while (i<len(time_array)-1):
    if lat_array[i] == lat_array[i+1] and long_array[i] == long_array[i+1]:
        bearing.append(0)
        i = i+1
    else:
        y = np.sin(difference(long_array)) * np.cos(lat_array[i+1])
        x = (np.cos(lat_array[i]) * np.sin(lat_array[i+1])) - (np.sin(lat_array[i]) * np.cos(lat_array[i+1]) * np.cos(long_array[i+1] - long_array[i]))
        theta = np.arctan2(y,x)
        if theta < 0:
            theta = theta + 2*np.pi
        bearing.append(theta)
        i = i+1

#Gradient of the bearing of the glider
rate = np.gradient(bearing)

#Plot bearing against time then select a section to zoom in on
"""zoom = [0, 0]

counter = 0
while counter < 2:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time_array, bearing)

    for i in range(0,1):
        cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

    x = coords[0]

    i = find_index(time_array, x)

    zoom[counter] = i

    counter = counter + 1

plt.plot(time_array[zoom[0]:zoom[1]], bearing[zoom[0]:zoom[1]])
plt.show()

plt.plot(time_array, rate)
plt.show()"""

dist_xyz = [0]
dist_xy = [0]
delta_z = [0]
R = 6371e3
i = 0
#Calculate distance as the crow flies in 3d space
while (i<len(time_array)-1):
    a = np.sin((difference(lat_array))/2)**2 + np.cos(lat_array[i]) * np.cos(lat_array[i+1]) * np.sin((difference(long_array))/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    d = R * c
    dz = difference(alt_gps_array)
    e = np.sqrt(d**2 +dz**2)
    delta_z.append(dz)
    dist_xy.append(d)
    dist_xyz.append(e)
    i = i+1

rel_lat = [x - lat_array[0] for x in lat_array]
x_pos = [(x*2*R)/(2*np.pi) for x in rel_lat]

rel_long = [y - long_array[0] for y in long_array]
y_pos = [(y*2*R)/(2*np.pi) for y in rel_long]

gps_minus_pres = [0]
i = 0
#Calculate difference in altitude between gps and pressure readings
while (i<len(time_array)-1):
    alt_dif = alt_gps_array[i] - alt_pres_array[i]
    gps_minus_pres.append(alt_dif)
    i = i+1
   
delta_t = [0]
i = 0
#Calculate the difference in time between each point
while (i<len(time_array)-1):
    dt = difference(time_array)
    delta_t.append(dt)
    i = i+1
   
ground_speed = [0]
i = 1
#Calculate the gorund speed as the crow flies
while (i<len(time_array)-1):
    speed = dist_xy[i]/delta_t[i]
    ground_speed.append(speed)
    i = i+1

height_speed = [0]
c = 1
#Calculate the rising/falling speed
"""while (c<len(time_array)-1): 
    speed = delta_z[c]/delta_t[c]
    height_speed.append(speed)
    c = c+1"""

"""plt.plot(y_pos, x_pos)
plt.show()"""

def thermaling():
    #Select the point which you want to check is in a thermal
    thermal = 0
    ia = 0
    while ia < 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Time (seconds since midnight)')
        ax.set_ylabel('Gps altitude (m)')
        ax.plot(time_array, alt_gps_array)

        for ib in range(0,1):
            cid = fig.canvas.mpl_connect('button_press_event', onclick)

        plt.show()

        x = coords[0]

        t = find_index(time_array, x)

        thermal = t
        
        ia = ia + 1

    i = 0
    while i<20:
        if np.abs(bearing[thermal] - bearing[thermal - 10+i]) > ((1)*np.pi):
            if thermal > thermal - 10+i:
                spiral = [thermal - 20, thermal + 20]
            else:
                spiral = [thermal - 10, thermal + 10]
        
            x_values = []
            y_values = []
            z_values = []
            i = spiral[0]
            #Plot the values for a specific thermal
            while (i<spiral[1]):
                j = ground_speed[i]*np.cos(bearing[i])
                k = ground_speed[i]*np.sin(bearing[i])
                l = alt_gps_array[i]
            
                x_values.append(j)
                y_values.append(k)
                z_values.append(l)
                i = i+1
            
            ax = plt.axes(projection = '3d')

            ax.scatter3D(x_values, y_values, z_values)
            ax.set_xlabel('East-West ground speed (m/s)')
            ax.set_ylabel('North-South ground speed (m/s)')
            ax.set_zlabel('Gps altitude (m)')
            plt.show()


            #Plot the actual thermal section of the flight in 3d
            rel_lat = [x - lat_array[0] for x in lat_array]
            x_pos = [(x*2*R)/(2*np.pi) for x in rel_lat]

            rel_long = [y - long_array[0] for y in long_array]
            y_pos = [(y*2*R)/(2*np.pi) for y in rel_long]

            ax = plt.axes(projection = '3d')

            ax.scatter3D(x_pos[spiral[0]:spiral[1]], y_pos[spiral[0]:spiral[1]], alt_gps_array[spiral[0]:spiral[1]])
            ax.set_xlabel('East-West')
            ax.set_ylabel('North-South')
            ax.set_zlabel('Gps altitude (m)')
            plt.show()

            
            #Plot the bearing against time of the thermal
            plt.plot(time_array[spiral[0]:spiral[1]], bearing[spiral[0]:spiral[1]])
            plt.xlabel('Time (seconds since midnight)')
            plt.ylabel('Bearing (radians)')
            plt.show()

            p = input("If you want to calculate the wind vector for this section input yes\n")

            if p == 'yes':

                #Using the Taubin-Newton method calculate the wind vector within the thermal
                wind_vector = Taubin_Newton(x_values, y_values)
                
                #Distance in x and y from each point to the center
                wind_vector_xrad = [abs(x - wind_vector[0]) for x in x_values]
                wind_vector_yrad = [abs(y - wind_vector[1]) for y in y_values]
                
                wind_vector_avg_xrad = np.average(wind_vector_xrad) #Average of the disatnces to the center
                wind_vector_avg_yrad = np.average(wind_vector_yrad)
                
                i = 0
                airspeed = []
                while i<len(x_values):
                    q = np.sqrt(wind_vector_xrad[i]**2 + wind_vector_yrad[i]**2)
                    airspeed.append(q)
                    i = i+1

                wind_vector_err = np.std(airspeed)/(np.sqrt(len(x_values)))
                knots_err = wind_vector_err*1.9438

                z = (alt_gps_array[spiral[1]] - alt_gps_array[spiral[0]])/(time_array[spiral[1]] - time_array[spiral[0]]) #Finding the average altitude velocity over the time selected
                wind_vector[2] = z
                wind_knots = wind_vector[3]*1.9438 #Calculating the wind speed in knots using a number found on Google
                
                wind_direction_rad = np.arctan2(wind_vector[1],wind_vector[0]) #Here the x and y values have been swapped to give the bearing relative to north
                #The % is used to return the remainder of the angle once divided by 2pi to give the angle as a positive value
                wind_direction_deg = rad_deg(wind_direction_rad) #Converting the wind direction to degrees to ease understanding

                print('wind vector is: [%5.2f, %5.2f, %5.2f]'% (wind_vector[0], wind_vector[1], wind_vector[2]))
                if(wind_direction_deg < 0):
                    wind_direction_deg = wind_direction_deg * (-1)
                    print('wind speed is: %5.2f +/-%5.2f m/s or %5.2f +/-%5.2f knots pointing %5.2f West of North'% (wind_vector[3], wind_vector_err, wind_knots, knots_err, wind_direction_deg))
                else:
                    print('wind speed is: %5.2f +/-%5.2f m/s or %5.2f +/-%5.2f knots pointing %5.2f East of North'% (wind_vector[3], wind_vector_err, wind_knots, knots_err, wind_direction_deg)) 

                wind_direction_deg = wind_direction_deg * (-1)

                #Plotting the wind vector
                ax = plt.axes(projection = '3d')
                ax.scatter3D(x_values, y_values, z_values)

                ax.set_xlabel('East-West')
                ax.set_ylabel('North-South')
                ax.set_zlabel('Altitude')
                ax.quiver(0, 0, np.average(z_values), wind_vector[0], wind_vector[1], wind_vector[2], color='r', length=wind_vector[3])
                plt.show()
                
                #Plot the actual thermal with the wind vector
                ax = plt.axes(projection = '3d')

                ax.scatter3D(x_pos[spiral[0]:spiral[1]], y_pos[spiral[0]:spiral[1]], alt_gps_array[spiral[0]:spiral[1]])
                
                ax.set_xlabel('East-West (m from start)')
                ax.set_ylabel('North-South (m from start)')
                ax.set_zlabel('Altitude (m)')
                ax.quiver(np.average(x_pos[spiral[0]:spiral[1]]), np.average(y_pos[spiral[0]:spiral[1]]), np.average(alt_gps_array[spiral[0]:spiral[1]]), wind_vector[0], wind_vector[1], wind_vector[2], color='r', length=wind_vector[3])
                plt.show()
                
                #Ploting the altitude time series with annotations
                ia = 0
                while ia < 1:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.set_xlabel('Time (seconds since midnight)')
                    ax.set_ylabel('Gps altitude (m)')
                    ax.plot(time_array, alt_gps_array)

                    for ib in range(0,1):
                        cid = fig.canvas.mpl_connect('button_press_event', onclick)

                    plt.show()

                    show = coords
                    
                    ia = ia + 1

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_xlabel('Time (seconds since midnight)')
                ax.set_ylabel('Gps altitude (m)')
                ax.plot(time_array, alt_gps_array)
                ax.plot(show[0], show[1], 'x')
                if (wind_direction_deg < 0):
                    wind_direction_deg = wind_direction_deg * (-1)
                    ax.annotate('[%5.2f, %5.2f, %5.2f]\n%5.2f +/-%5.2f m/s or %5.2f +/-%5.2f knots\nPointing %5.2f° West of North'% (wind_vector[0], wind_vector[1], wind_vector[2], wind_vector[3], wind_vector_err, wind_knots, knots_err, wind_direction_deg), xy=(show[0], show[1]), xytext=(show[0]+100, show[1]-300), fontsize=12)
                else:
                    ax.annotate('[%5.2f, %5.2f, %5.2f]\n%5.2f +/-%5.2f m/s or %5.2f +/-%5.2f knots\nPointing %5.2f° East of North'% (wind_vector[0], wind_vector[1], wind_vector[2], wind_vector[3], wind_vector_err, wind_knots, knots_err, wind_direction_deg), xy=(show[0], show[1]), xytext=(show[0]+100, show[1]-300), fontsize=12)
                plt.show()

                #Colour application
                my_colour = np.where(np.gradient(alt_gps_array)>=0, 'red', 'blue')

                """ax = plt.axes(projection = '3d')

                ax.scatter3D(x_pos, y_pos, alt_gps_array, color = my_colour)
                plt.show()"""
                break
            else:
                go = thermaling()
        else:
            if i == 19:
                print('The thermal method does not work here, try again')
                go = thermaling()
            else:
                i = i+1

    return 0

go = thermaling()