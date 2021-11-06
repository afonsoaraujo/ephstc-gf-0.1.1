import os, sys
os.chdir('/home/vitor/Desktop/MOD16')
from ftplib import FTP, error_perm

def walk_dir(f, dirpath):
    original_dir = f.pwd()
    try:
        f.cwd(dirpath)
    except error_perm:
        return # ignore non-directories and ones we cannot enter
    names = f.nlst()
    names.sort
    for i in range(len(names)):
        f.cwd(names[i]) # return to cwd of our caller
        #######################################################################################################
        arquivo1 = f.nlst('*h13v11*')
        print "Donwload em andamento do arquivo: ",arquivo1[0] 
        f.retrbinary("RETR " + arquivo1[0],open(arquivo1[0], 'wb').write)
        print "Download concluido"
        print
        print
        #######################################################################################################
        arquivo2 = f.nlst('*h14v11*')
        print "Donwload em andamento do arquivo: ",arquivo2[0] 
        f.retrbinary("RETR " + arquivo2[0],open(arquivo2[0], 'wb').write)
        print "Download concluido"
        print
        print
        f.cwd(dirpath)

print "##################################################################"
print "##################################################################"
print "######                                                      ######"
print "######               DOWNLOAD DE DADOS MODIS                ######"
print "######                       MOD16A2                        ######"
print "######                                                      ######"
print "######                   Desenvolvido por                   ######"
print "######              Lab. Hidrologia-COOPE-UFRJ              ######"
print "######                 site:ftp.ntsg.umt.edu                ######"
print "######                                                      ######"
print "##################################################################"
print "##################################################################"

print
print
print
print "##################################################################"
print "INICIALIZANDO O PROGRAMA"
print "##################################################################"
print
print
print

f = FTP('ftp.ntsg.umt.edu')
f.login()
walk_dir(f,'/pub/MODIS/NTSG_Products/MOD16/MOD16A2.105_MERRAGMAO/Y2000')
walk_dir(f,'/pub/MODIS/NTSG_Products/MOD16/MOD16A2.105_MERRAGMAO/Y2001')
walk_dir(f,'/pub/MODIS/NTSG_Products/MOD16/MOD16A2.105_MERRAGMAO/Y2002')
walk_dir(f,'/pub/MODIS/NTSG_Products/MOD16/MOD16A2.105_MERRAGMAO/Y2003')
walk_dir(f,'/pub/MODIS/NTSG_Products/MOD16/MOD16A2.105_MERRAGMAO/Y2004')

print
print
print
print "##################################################################"
print "O PROGRAMA TERMINOU COM SUCESSO"
print "##################################################################"
print
print
print
