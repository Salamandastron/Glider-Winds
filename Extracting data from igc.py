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