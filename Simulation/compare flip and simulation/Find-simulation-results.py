import csv, math, time


num_files = 49 #number of files to be scanned
dt = 0.2;           #Step size
t_frame = 1.6      #time between each frame
t_bleach = 0.4
t_to_bleach = 0     # t_to_bleach * t_frame = first bleach
dir = "../Mesh/solutions_logistic/"
nucleus = []
cyto = []
nucleus1 = []
cyto1 = []
indices = []
nucleus2 = []
cyto2 = []
indices2 = []


for i in range(0,num_files):
    if i >= t_frame/dt*t_to_bleach-2:
        # Read from file
        if i >= 1000:
            searchfile = open(dir + "c00" + str(i) + ".vtu", "r")
        elif i >= 100:
            searchfile = open(dir + "c000" + str(i) + ".vtu", "r")
        elif i >= 10:
            searchfile = open(dir + "c0000" + str(i) + ".vtu", "r")
        else:
            searchfile = open(dir + "c00000" + str(i) + ".vtu", "r")
                
        for line in searchfile:
            if "Name=\"u\"" in line:
                break
        searchfile.close()

        # Just before bleaching
        if  t_frame/dt*math.ceil((i+1)*dt/t_frame)-2 == i:
            print t_frame/dt*math.ceil((i+1)*dt/t_frame)-2
                    
            nums = line[53:-14].split( )
            nums = map(float, nums)
            nucleus.append(max(nums))
            cyto.append(min(nums))
        #        nucleus.append(sum(nums[0:250])/len(nums[0:250]))
        #        cyto.append(sum(nums[len(nums)-250:len(nums)])/len(nums[0:250]))
        # Straight after bleching
        elif t_frame/dt*math.ceil((i-3)*dt/t_frame)+t_bleach/dt-1 == i:
                
            nums1 = line[53:-14].split( )
            nums1 = map(float, nums1)
            nucleus1.append(max(nums1))
            cyto1.append(min(nums1))
        #        nucleus1.append(sum(nums1[0:250])/len(nums1[0:250]))
        #        cyto1.append(sum(nums1[len(nums1)-250:len(nums1)])/len(nums1[0:250]))
            indices.append(i)
        # In between two bleaches
        elif i == int((math.floor(i/(t_frame/dt))*t_frame + t_bleach + math.floor(t_frame*0.5))/dt+1):

            nums2 = line[53:-14].split( )
            nums2 = map(float, nums2)
            nucleus2.append(max(nums2))
            cyto2.append(min(nums2))
        #        nucleus2.append(sum(nums2[0:250])/len(nums2[0:250]))
        #        cyto2.append(sum(nums2[len(nums2)-250:len(nums2)])/len(nums2[0:250]))
            indices2.append(i)


for item in indices:
    print item
for item in indices2:
    print item
# open a file for writing.
csv_out = open('simulation-results-before-bleach.csv', 'wb')
mywriter = csv.writer(csv_out)

# writerow - one row of data at a time.
for row in zip(cyto, nucleus):
    mywriter.writerow(row)

csv_out.close()


# open a file for writing.
csv_out = open('simulation-results-after-bleach.csv', 'wb')
mywriter = csv.writer(csv_out)

# writerow - one row of data at a time.
for row in zip(cyto1, nucleus1):
    mywriter.writerow(row)

csv_out.close()


# open a file for writing.
csv_out = open('simulation-results-between-bleaches.csv', 'wb')
mywriter = csv.writer(csv_out)

# writerow - one row of data at a time.
for row in zip(cyto2, nucleus2):
    mywriter.writerow(row)

csv_out.close()
