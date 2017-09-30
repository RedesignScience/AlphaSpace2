"""
Usage:
put this in the folder where you have alphaspace as sub-directories, then run the program

for instance:
python change_lig_resid.py -8


"""




import os,sys



def space_fill(input,length=0):
    output = ""
    string = str(input)
    if length != 0 and length > len(string):
        for i in range(length-len(string)):
            output+= " "
        output+= string
    else:
        output = string
    return output




def update_lig_resid(directories,adjustment = -8):
    for folder in directories:
        as_folders = [os.path.join(folder, d) for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]
        for as_folder in as_folders:
            lig_path = as_folder+"/pdb_out/lig.pdb"
            try:
                with open(lig_path,'r') as handle:
                    lines = handle.readlines()
                new_lines = []
                for line in lines:
                    if len(line.split()) > 1:
                        print(line)
                        new_lines.append(line[:8]+ space_fill(int(line[8:11])+adjustment,length=3)+line[11:])
                        print(new_lines[-1])
                with open(lig_path,'w') as handle:
                    handle.writelines(lines)
            except:
                pass








if __name__ == '__main__':
    try:
        adjustment = int(sys.argv[1])
    except:
        raise "please specify adjustment value"
    cwd = os.getcwd()
    directories = [os.path.join(cwd, d) for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]
    update_lig_resid(directories,adjustment)