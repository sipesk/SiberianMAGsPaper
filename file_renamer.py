import os
# uses os module https://docs.python.org/2/library/os.html
top_dir = os.getcwd()
while True:
    mydir = 'all_fnas'
    try:
        os.makedirs(mydir)
        break
    except Exception as e:
        if e.errno != os.errno.EEXIST:
            raise   
        # time.sleep might help here
        pass
print('TopDir =')
print(top_dir)

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def get_immediate_subfiles(a_dir):
    return [name for name in os.listdir(a_dir)
            if not os.path.isdir(os.path.join(a_dir, name))]

second_level_dirs = get_immediate_subdirectories(top_dir)

for directory in second_level_dirs:
    print(get_immediate_subfiles(directory))
    subfiles = get_immediate_subfiles(directory)

    for subfileX in subfiles:
        name_without_ext = ''
        # if subfileX.endswith('.txt'):
        #     print('instance of txt file')
        #     name_without_ext = subfileX.replace('.txt', '', 1)
        #     print(name_without_ext)
        # elif subfileX.endswith('.text'):
        #     print('instance of text file')
        #     name_without_ext = subfileX.replace('.text', '', 1)
        #     print(name_without_ext)
        if subfileX.endswith('.fna'):
            if not subfileX.endswith('.genes.fna') and not subfileX.endswith('.intergenic.fna'):
                print(directory)
                print('instance of fna file')
                print(subfileX)
                file_obj_fna = open(top_dir + '/' + directory + '/' + subfileX, 'r')
                first_line = file_obj_fna.readline()
                first_line = first_line.split(':')[0]
                file_obj_fna.close()
                new_name = first_line.replace('>', '', 1)
                new_name = new_name.replace(' ', '_')
                os.rename(directory, new_name)
                # os.rename(str(top_dir) + '/' + str(directory) + '/' + str(subfileX), str(top_dir) + '/all_fnas/' + str(new_name))

# end of for loop
all_fnas_dir = top_dir + '/' + 'all_fnas'
all_fna_files = get_immediate_subfiles(all_fnas_dir)

#for fna in all_fna_files:
#     os.replace(top_dir + '/' + 'all_fnas/' + fna, top_dir + '/' + 'all_fnas/' + fna + '.fa')

# for fna in all_fna_files:
#     print(fna)
#     file_obj = open(top_dir + '/' + 'all_fnas/' + str(fna), 'r')
#     text = file_obj.read()
#     file_obj.close()
#     file_obj = open(top_dir + '/' + 'all_fnas/' + str(fna), 'w')
#     myStr = text.encode('ascii',errors='ignore').decode()
#     file_obj.write(myStr)
#     file_obj.close()
