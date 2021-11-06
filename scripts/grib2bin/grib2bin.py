#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

##############################################################
# Realiza download do Conjunto de Dados GLDAS_NOAH025SUBP_3H
##############################################################

#BIBLIOTECAS
from sys import argv, exit, stdout
from subprocess import Popen,PIPE
from os import path,remove,rename
from datetime import datetime, date, time, timedelta
from copy import copy
from pydap.client import open_url
import sqlite3
sqlite = sqlite3

try:
 import pygtk
 pygtk.require("2.00")
except:
 pass
try:
 import gtk
except:
 print("GTK Not Availible")
 exit(1)

#----------------
#GLOBAL VARIABLES
#----------------
url = "http://agdisc.gsfc.nasa.gov/dods/GLDAS_NOAH025SUBP_3H"
fname_db = 'config.db'

#-----------------
# CTL INFORMATION 
#-----------------

class ctl_info:
 def __init__( self ):
  self.is_empty = True
  #General Information
  self.filename = str()
  self.title = str()
  self.undef = float()
  self.options = str()
 
  #Spatial Information
  #Longitude
  self.nx = int()
  self.xtype = str()
  self.longitude = list()
  self.dx = float()
 
  #Latitude
  self.ny = int()
  self.ytype = str()
  self.latitude = list()
  self.dy = float()
 
  #Level
  self.nz = int()
  self.ytype = str()
  self.levels = list()
  self.dz = float()
 
  #Time
  self.nt = int()
  self.ttype = str()
  self.times = list()
  self.dt = str()
  self.udt = str()
 
  #Variables
  self.nvar = int()
  self.vars_name = list()
  self.var_description = list()

 def fromfile( self, filename = str() ):
  f = open(filename,'rt')
  ctl = f.readlines()
  f.close()
  if ctl[-1] == '\n': del(ctl[-1])
  #Read CTL File
  line_is_var =  False
  for n,line in enumerate(ctl):
   words = line.split()
   if words[0].lower() == 'dset':
    self.filename = words[1]
   if words[0].lower() == 'title':
    for word in words[1:-1]:
     self.title += word + ' '
    self.title += words[-1]
   elif words[0].lower() == 'undef':
    self.undef = float(words[1])
   elif words[0].lower() == 'options':
    self.options = line[6:]
   elif words[0].lower() == 'xdef':
    self.nx = int(words[1])
    self.xtype = words[2]
    self.dx = float(words[4])
    if self.xtype.lower() == 'linear':
     self.longitude.append( float(words[3]) )
     i = 1
     while i != self.nx:
      self.longitude.append( self.longitude[0] + i * self.dx )
      i += 1
   elif words[0].lower() == 'ydef':
    self.ny = int(words[1])
    self.ytype = words[2]
    self.first_lat = float(words[3])
    self.dy = float(words[4])
    if self.ytype.lower() == 'linear':
     self.latitude.append( float(words[3]) )
     i = 1
     while i != self.ny:
      self.latitude.append( self.latitude[0] + i * self.dy )
      i += 1
   elif words[0].lower() == 'zdef':
    self.nz = int(words[1])
    self.ztype = words[2]
    self.levels.append(float(words[3]))
    self.dz = float(words[4])
   elif words[0].lower() == 'tdef':
    self.nt = int(words[1])
    self.ttype = words[2]
    self.first_time = datetime.strptime(words[3].lower(), "%Hz%d%b%Y")
    for character in words[4]:
     if character.isdigit():
      self.dt = self.dt + character
     else:
      self.udt = self.udt + character #FIX-ME VALIDAR DADOS
    if self.udt.lower() == 'mn': self.dt = timedelta(minutes = int(self.dt))
    if self.udt.lower() == 'hr': self.dt = timedelta(hours = int(self.dt))
    if self.udt.lower() == 'dy': self.dt = timedelta(days = int(self.dt))
    #IMPLEMENTAR MESES
    
   elif words[0].lower() == 'vars':
    self.nvar = int(words[1])
    line_is_var =  True
    continue
   elif words[0].lower() == 'endvars':
    line_is_var =  False
   if line_is_var:
    self.vars_name.append( words[0].split('=>')[0] )
    self.var_description.append( line.split('**')[-1] )
  self.is_empty = False
 def tofile( self, fname = str() ):
  d = ' '
  if self.is_empty:
   print 'Error writing Ctl File: Ctl information is empty'
  else:
   if fname[(len(fname)-3):].lower() != 'ctl': fname = fname + '.ctl'
   f = open(fname,'w')
   f.write('dset '+ self.filename + '\n')   
   f.write('options '+ self.options + '\n')
   f.write('title '+ self.title + '\n')
   f.write('undef '+ str(self.undef) + '\n')
   if self.xtype.lower() == 'linear':
    f.write('xdef '+ str(self.nx) + d + self.xtype + d + str(self.longitude[0]) + d + str(self.dx)+'\n')
   if self.ytype.lower() == 'linear':
    f.write('ydef '+ str(self.ny) + d + self.ytype + d + str(self.latitude[0]) + d + str(self.dy)+'\n')
   if self.ztype.lower() == 'linear':
    f.write('zdef '+ str(self.nz) + d + self.ztype + d + str(self.levels[0]) + d + str(self.dz)+'\n')
   if self.ttype.lower() == 'linear':
    if self.udt.lower() == 'mn':
     f.write('tdef '+ str(self.nt) + d + self.ttype + d + self.first_time.strftime("%Hz%d%b%Y") + d + str(int(self.dt.seconds/60.0))+self.udt+'\n')
    if self.udt.lower() == 'hr':
     f.write('tdef '+ str(self.nt) + d + self.ttype + d + self.first_time.strftime("%Hz%d%b%Y") + d + str(int(self.dt.seconds/3600.0))+self.udt+'\n')
    if self.udt.lower() == 'dy':
     f.write('tdef '+ str(self.nt) + d + self.ttype + d + self.first_time.strftime("%Hz%d%b%Y") + d + str(int(self.dt.days))+self.udt+'\n')
   f.write('vars ' + str(self.nvar) + '\n')
   for nd, descript in enumerate(self.var_description):
    f.write(d + self.vars_name[nd] + ' 0 99 '+ descript)
   f.write('\nendvars\n')
   f.close()
   
 def fromurl(self, url = str()):
  dataset = open_url(url)
  self.dataset = dataset
  self.filename = copy(url)
  self.title = copy(dataset.attributes['NC_GLOBAL']['title'])
  self.undef = copy(dataset[dataset.keys()[0]].attributes['missing_value'])
  self.nt,self.ny,self.nx = copy(dataset[dataset.keys()[0]].shape)

  self.xtype = 'linear'
  self.ytype = 'linear'
  self.ztype = 'linear'
  self.ttype = 'linear'

  longitude = dataset.lon.attributes['minimum']
  self.dx = int((dataset.lon.attributes['maximum']-longitude+1)/self.nx*100)/100.0
  for nlon in range(self.nx): self.longitude.append(longitude+nlon*self.dx)

  latitude = dataset.lat.attributes['minimum']
  self.dy = int((dataset.lat.attributes['maximum']-latitude+1)/self.ny*100)/100.0
  for nlat in range(self.ny): self.latitude.append(latitude+nlat*self.dy)

  self.dz = 1
  self.levels.append(1)
 
  self.first_time = datetime.strptime(dataset.time.attributes['minimum'].lower(), "%Hz%d%b%Y")
  self.last_time = datetime.strptime(dataset.time.attributes['maximum'].lower(), "%Hz%d%b%Y")
  #NEED TO FIX - PROBLEM: NEED TO BE MORE GENERAL
  self.dt = (self.last_time - self.first_time)/self.nt
  self.dt = timedelta( hours = round(self.dt.seconds/3600.0) )
  self.udt = 'hr'

  self.nvar = len(dataset.keys())-3
  self.vars_name = copy(dataset.keys()[0:-3])
  for nvar in range(self.nvar):
   self.var_description.append( dataset[dataset.keys()[nvar]].attributes['long_name'][3:] )
  self.is_empty = False
#------------------
#GTK WINDOWS CLASS
#------------------

class main_Window:
 def __init__( self ):
  self.builder = gtk.Builder()
  self.builder.add_from_file( "main.glade" )

  self.window = self.builder.get_object('mainWindow')
  self.window.connect('destroy', lambda w: exit(0))
  self.txt_input = self.builder.get_object("txt_input")
  self.txt_input.set_text(url)
  self.txt_prefix = self.builder.get_object("txt_prefix")
#Connection with database
  if path.exists(fname_db):
   self.connection = sqlite.connect(fname_db)
   self.cursor = self.connection.cursor()
  else:
   self.connection = sqlite.connect(fname_db)
   self.cursor = self.connection.cursor()
   self.cursor.execute('CREATE TABLE configs ( \
   id INTEGER PRIMARY KEY, \
   name TEXT, \
   prefix TEXT, \
   firstlon REAL, \
   lastlon REAL, \
   firstlat REAL, \
   lastlat REAL, \
   firsttime INTEGER, \
   lasttime INTEGER )')
   self.cursor.execute("""INSERT INTO configs VALUES (\
   null, 'South America','sa',-83.75,-31.25,-53.75,16.25,1,1)""")
   self.cursor.execute("""INSERT INTO configs VALUES (\
   null, 'Brazil','br',-75.00,-34.00,-35.00,5.00,1,1)""")
   self.cursor.execute("""INSERT INTO configs VALUES (\
   null, 'Negro River','nr',-74.00,-58.30,-3.30,5.48,1,1)""")
   self.connection.commit()

  self.cbox_config = self.builder.get_object("cbox_config")
  self.liststore_config = gtk.ListStore(str)
  cell = gtk.CellRendererText()
  self.cbox_config.pack_start(cell)
  self.cbox_config.add_attribute(cell, 'text', 0)
  self.liststore_config.append(['Server Information'])
  self.cursor.execute('SELECT * FROM configs')
  self.list_of_configs = self.cursor.fetchall()
  
  for item in self.list_of_configs:
   self.liststore_config.append([item[1]])

  self.cbox_config.set_model(self.liststore_config)
  self.cbox_config.set_active(0)
  self.cbox_config.connect( "changed", lambda w: self.load_config() )
   
  self.btt_new = self.builder.get_object("btt_new")
  self.btt_new.connect( "clicked", lambda w: self.new_config() )
  self.btt_save = self.builder.get_object("btt_save")
  self.btt_save.connect( "clicked", lambda w: self.save_config() )
  self.btt_delete = self.builder.get_object("btt_delete")
  self.btt_delete.connect("clicked", lambda w: self.remove_config() )

  self.first_lat = self.builder.get_object("spin_first_lat")
  self.last_lat = self.builder.get_object("spin_last_lat")
  self.first_lon = self.builder.get_object("spin_first_lon")
  self.last_lon = self.builder.get_object("spin_last_lon")
  self.spin_first_time = self.builder.get_object("spin_first_time")
  self.spin_first_time.connect("value-changed", lambda w: self.update_date() )
  self.spin_last_time = self.builder.get_object("spin_last_time")
  self.spin_last_time.connect("value-changed", lambda w: self.update_date() )

  self.cbox_var = self.builder.get_object("cbox_var")
  self.progressbar = self.builder.get_object("progressbar")

  self.btt_quit = self.builder.get_object("btt_quit")
  self.btt_quit.connect("clicked", lambda w: exit(0))
  self.btt_download = self.builder.get_object("btt_download")
  self.btt_download.connect("released", lambda w: self.download() )

  self.info = ctl_info()
  self.load_config_from_ctl_info()
  self.info.fromurl(url)#IMPLEMENT IDENTIFY IF IS A URL
  self.update_date()
  gtk.main()
  
 def update_date( self ):
  self.lbl_fist_time = self.builder.get_object("lbl_first_time")
  self.lbl_fist_time.set_label( \
  (self.info.first_time + (int(self.spin_first_time.get_value())-1)*self.info.dt).strftime("%HZ%d%b%Y") )

  self.lbl_last_time = self.builder.get_object("lbl_last_time")
  self.lbl_last_time.set_label(\
  (self.info.first_time + (int(self.spin_last_time.get_value())-1)*self.info.dt).strftime("%HZ%d%b%Y") )

 def download( self ):
  self.progressbar.set_fraction(0.0)
  lon_choosed = str( self.first_lon.get_value() ) + " " + str( self.last_lon.get_value() )
  lat_choosed = str( self.first_lat.get_value() ) + " " + str( self.last_lat.get_value() )
  var_name = self.info.vars_name[ self.cbox_var.get_active() ]
  first_time = int( self.spin_first_time.get_value() )
  last_time = int( self.spin_last_time.get_value() )

  #UPDATE THE CTL INFORMATION
  #NEED FIXES: 1) consider cases of non linear x,y,z,t levels

  var_info = copy(self.info)
  var_info.filename = "^"+self.txt_prefix.get_text()+"_" + var_name +'_%y4%m2%d2%h2.r4'
  var_info.options = 'template'
  var_info.nx = int( (self.last_lon.get_value() - self.first_lon.get_value() ) / var_info.dx ) + 1
  var_info.longitude = list()
  for i in range(var_info.nx):
   var_info.longitude.append(self.first_lon.get_value()+i*var_info.dx)

  var_info.ny = int( (self.last_lat.get_value() - self.first_lat.get_value() ) / var_info.dy ) + 1
  var_info.latitude = list()
  for i in range(var_info.ny):
   var_info.latitude.append(self.first_lat.get_value()+i*var_info.dy)

  var_info.nt = int( self.spin_last_time.get_value() - self.spin_first_time.get_value() ) + 1
  var_info.first_time = self.info.first_time + int(self.spin_first_time.get_value() - 1) * var_info.dt
  
  var_info.nvar = 1
  var_info.vars_name = list()
  var_info.var_description = list()
  var_info.vars_name.append(self.info.vars_name[ self.cbox_var.get_active() ])
  var_info.var_description.append( self.info.var_description[ self.cbox_var.get_active() ])
  var_info.tofile(self.txt_prefix.get_text()+'_' + var_name + '.ctl')
  var = var_info.dataset[var_info.vars_name[0]]
#WRITE BINARY file

  nflat = self.info.latitude.index(var_info.latitude[0])
  nllat = self.info.latitude.index(var_info.latitude[-1]) + 1
  nflon = self.info.longitude.index(var_info.longitude[0])
  nllon = self.info.longitude.index(var_info.longitude[-1]) + 1
  for time_choosed in range( first_time, last_time + 1 ):
   print time_choosed
   current_date = self.info.first_time + int(time_choosed - 1) * var_info.dt
   fname = self.txt_prefix.get_text()+'_'+var_name+"_"+current_date.strftime('%y%m%d%H')+'.r4'
   matrix = var.array[time_choosed-1,nflat:nllat,nflon:nllon].astype('=f4')
   matrix.tofile(fname)
   self.progressbar.set_fraction(float(time_choosed)/var_info.nt)
  print 'Download complete.'
 def save_config( self ):
  configs = self.cbox_config.get_active()
  print configs
  #self.cursor.execute('INSERT INTO configs VALUES (null,' + configs)
  #self.connection.commit()
 def load_config_from_ctl_info( self ):
  self.info.fromurl(url)
  #LONGITUDE
  adjustment = gtk.Adjustment(self.info.longitude[0],self.info.longitude[0],self.info.longitude[-1], self.info.dy)
  self.first_lon.set_adjustment( adjustment )
  self.first_lon.set_value( self.info.longitude[0] )

  adjustment = gtk.Adjustment(self.info.longitude[0],self.info.longitude[0],self.info.longitude[-1], self.info.dy)
  self.last_lon.set_adjustment( adjustment )
  self.last_lon.set_value( self.info.longitude[-1] )

  #LATITUDE
  adjustment = gtk.Adjustment(self.info.latitude[0],self.info.latitude[0],self.info.latitude[-1], self.info.dx)
  self.first_lat.set_adjustment( adjustment )
  self.first_lat.set_value( self.info.latitude[0] )
  adjustment = gtk.Adjustment(self.info.latitude[0],self.info.latitude[0],self.info.latitude[-1], self.info.dx)
  self.last_lat.set_adjustment( adjustment )
  self.last_lat.set_value( self.info.latitude[-1] )

  #TIME
  adjustment = gtk.Adjustment( 1, 1, self.info.nt, 1 )
  self.spin_first_time.set_adjustment( adjustment )
  self.spin_first_time.set_value( 1 )

  adjustment = gtk.Adjustment( self.info.nt, 1, self.info.nt, 1 )
  self.spin_last_time.set_adjustment( adjustment )
  self.spin_last_time.set_value( self.info.nt )

  #VAR
  liststore = gtk.ListStore(str)
  cell = gtk.CellRendererText()
  self.cbox_var.pack_start(cell)
  self.cbox_var.add_attribute(cell, 'text', 0)
  for titulo in self.info.var_description: liststore.append([titulo])
  self.cbox_var.set_model(liststore)
  self.cbox_var.set_active(0)

 def load_config( self ):
  self.cursor.execute('SELECT * FROM configs')
  self.list_of_configs = self.cursor.fetchall()
  n = self.cbox_config.get_active()
  if n == 0 :
   pass
  else:
   item = self.list_of_configs[n-1]

   self.txt_prefix.set_text( item[2] )
   self.first_lon.set_value( item[3] )
   self.last_lon.set_value( item[4] )   
   self.first_lat.set_value( item[5] )
   self.last_lat.set_value( item[6] )
   self.spin_first_time.set_value( item[7] )
   self.spin_last_time.set_value( item[8] )

   self.first_lon.update()
   self.last_lon.update()
   self.first_lat.update()
   self.last_lat.update()
   self.spin_first_time.update()
   self.spin_last_time.update()
   
 def new_config( self ):
  pass
 def remove_config( self ):
  pass

main_Window()
