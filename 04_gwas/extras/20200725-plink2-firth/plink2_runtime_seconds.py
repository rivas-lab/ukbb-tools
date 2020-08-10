import sys

def plink2_runtime_seconds(log_f):
    from datetime import datetime
    
    with open(log_f) as f:
        lines = f.read().splitlines()
        
    times = [
        datetime.strptime(
            x.replace('Start time: ', '').replace('End time: ', ''),
            '%c'
        ) for x in list(filter(
            lambda l: l.startswith('Start time') or l.startswith('End time'), 
            lines
        ))
    ]
    
    return(sum([ 
        int(( times[2*i+1] - times[2*i] ).total_seconds()) 
        for i in range(int(len(times) / 2))
    ]))


if __name__ == '__main__':
    print(plink2_runtime_seconds(str(sys.argv[1])))
