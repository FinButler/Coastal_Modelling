import csv
def val(t):
        InputWaveReader = csv.reader(open('raw_data/InputWave.csv', 'r'), delimiter='\t')
        data=[]
        for (time, height) in InputWaveReader:
                data.append((float(time), float(height)))
                
        for i in range(1,len(data)):
                if data[i][0]<t:
                        continue
                t1=data[max(0,i-1)][0]
                t2=data[i][0]
                h1=data[max(0,i-1)][1]
                h2=data[i][1]
                return h1*(t-t2)/(t1-t2)+h2*(t-t1)/(t2-t1)
        
        return data[-1][1]
