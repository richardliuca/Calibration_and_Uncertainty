import os, pyexcel
from pathlib import Path

def lvm_to_csv(workDir = None, allfiles = False, target_file = None,
            newName = 'Angle_of_Attack', ignore = None):
    ''' Converts LVM files without headers to csv files
        This function is optimized for special cases
    '''
    try:
        workDir = Path(workDir)
        os.chdir(workDir)
        filesList = [file for file in os.listdir(workDir) if '.lvm' in file]
        if not ignore == None and ignore in filesList:
            filesList.remove(ignore)
        elif not ignore == None and not ignore in filesList:
            print('Ignore file does not exist, not ignoring')
        else:
            pass
    except:
        return print('The work directory does not exist')

    if allfiles == False and target_file in filesList:
        filesList = [target_file]
    elif allfiles == False and not target_file in filesList:
        return print('Target file is not in directory')
    elif allfiles == True:
        print('Reading all files')
    else:
        return print('I got confused')

    for file in filesList:
        print('Reading file: {}'.format(file))
        # Opening LVM file and read the data as string
        f = open(file, 'r+')
        content = f.read()
        f.close()
        # Find the last instance of End of Header
        substring = 'End_of_Header'
        lastHeader = content.rfind(substring)
        # Find the last new line after the last header
        lastEnter = content.find('\n', lastHeader)
        # Take all string from the last new line after the laster header
        # until the end of file
        newContent = content[lastEnter+1:-1]
        # Separate entire file string by row (new line)
        contentSplit = newContent.split('\n')
        # Separate the column title from the first entries of the content
        #  and split by tab
        columnTitle = contentSplit[0].split('\t')
        # NOTE: This is made fo personal cases where the last column title
        # is removed and the second to last column title is renamed
        columnTitle.pop()
        columnTitle[-1] = newName
        # Preparaing the new list for properly sectioned file stream
        columnData = [columnTitle]
        # Iterate through the rows of the file contente
        for i in range(len(contentSplit)-1):
            # For each row split the content by tab indication starting
            # from the second row (frist row is column title)
            row = contentSplit[i+1].split('\t')
            # Remove last entries if it is larger than length of 7
            [row.pop() for _ in range((len(row)- 7))]
            if i == 0:
                # Look for the first entries of last column
                alpha = row[-1]
            elif '' in row and row[-1] == '':
                # Fill all empty entries and replace wth the first entries
                row[-1] = alpha
            else:
                pass
            # Building 2d list with each row
            columnData.append(row)
        # Save file as .csv with the same name
        newFileName = file.replace('.lvm', '.csv')
        pyexcel.save_as(array = columnData, dest_file_name = newFileName)

if __name__ == '__main__':
    dir = r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\10.26.2018\Angle Attack Cali'
    settings = {'workDir': dir, 'allfiles': True, 'ignore':'angleofattackcal_14.lvm'}
    lvm_to_csv(**settings)
