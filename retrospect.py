#!/usr/bin/python
from numpy import random, sqrt,log
import sys
import re
import MySQLdb
import numpy as np
from random import randrange, uniform
import copy
#mind to clear cache when using BM function on different column
class retro_cache:
  bm_num1 = None
  bm_num2 = None
  bm_prev_flg = 0
  table_name = ""

def clear_retro_cache():
  retro_cache.bm_num1 = None
  retro_cache.bm_num2 = None
  retro_cache.bm_prev_flg = 0
  retro_cache.table_name = ""


def boxmullerpolar():
    if (retro_cache.bm_prev_flg):
      retro_cache.bm_prev_flg = 0
      return retro_cache.bm_num2

    while True:
      # uniformly distributed values between 0 and 1
      xt = random.rand()
      yt = random.rand()
      x = (2 * xt) - 1
      y = (2 * yt) - 1
      r = x*x + y*y
      if (r< 1 and r > 0): break

    pl = sqrt(-2.0*log(r)/r)
    retro_cache.bm_num1 = pl * x;
    retro_cache.bm_num2 = pl * y;
    retro_cache.bm_prev_flg = 1
    return retro_cache.bm_num1

def gen_frm_equ(x,cf):
    poly = np.poly1d(cf)
    tmp = np.polyval(poly, x)
    return tmp

def find_polyfit(x_arr,y_arr): 
    return np.polyfit(x_arr, y_arr, deg=1)#you can change degree as you see fit

def trnsfrm_bmp_2imputed(bmp,mean,std):
    tmp = bmp*std + mean
    return tmp

def my_mean_squared_error(org,pred,len_org,len_pred):
  err = org - pred
  err = err * err
  #l_diff = len_org - len_pred ERROR 1
  l_diff = len_org - 2
  s_err = sum(err)
  if l_diff == 0:
    return s_err
  else:
    return s_err/(abs(l_diff))

def find_noise(perr):
  l = len(perr)
  i = 0
  noise = []
  while (i < l):
    noise.append(trnsfrm_bmp_2imputed( boxmullerpolar(),0,sqrt(perr[i])))
    i = i+1
  return noise

def strip_list(lst,missing_list):
  del lst[missing_list[0]-1:missing_list[-1]]
  return lst

def compute_ir(xa,mean):
  ir = 0
  for i in xrange(0,len(xa)):
    #ir = ir + (np.exp(xa[i])-mean)*(np.exp(xa[i])-mean) ERROR 2
    ir = ir + ((xa[i]-mean)*(xa[i]-mean))
  return ir

# transformation function
def find_missingvalues( idep, dep, missing_list):
  genl = []
  striped_idep = copy.deepcopy(idep)
  striped_dep = copy.deepcopy(dep)
  striped_idep= strip_list( striped_idep, missing_list)
  striped_dep = strip_list( striped_dep, missing_list)
  striped_idep.sort()
  striped_dep.sort()
  
  xa = np.array(striped_idep[:])
  xa = np.log(xa)
  ya = np.asarray(striped_dep[:])
  ya = np.log(ya)
  cf = np.polyfit(xa, ya, deg=1)#you can change degree as you see fit
  
  xaf = np.array(idep)
  xaf = np.log(xaf)
  yaf = gen_frm_equ(xaf,cf)
  pred_dep = list(yaf)
  
  striped_pred_dep = copy.deepcopy(pred_dep)
  striped_pred_dep = strip_list( striped_pred_dep, missing_list)
  ya_p_s = np.array(striped_pred_dep)
  mse = my_mean_squared_error(ya,ya_p_s,len(striped_dep),len(missing_list))
  d_m = np.mean(xa)
  #print xa, np.exp(xa)
   #=$V$2*(1+(1/3)+((A3-$X$2)*(A3-$X$2)/(($A$2-$X$2)*($A$2-$X$2)+($A$5-$X$2)*($A$5-$X$2)+($A$6-$X$2)*($A$6-$X$2))))
  p_err = []
  ir = compute_ir(xa,d_m)
  l = len(xa)
  for i in missing_list:
    #t = mse*(1+(1/l)+((idep[i-1]/1000-d_m)*(idep[i-1]/1000-d_m)/ir))
    #t = mse*(1+(1/l)+((idep[i-1]-d_m)*(idep[i-1]-d_m)/ir)) ERROR 3
    t = mse*(1+(1/l)+((log(idep[i-1])-d_m)*(log(idep[i-1])-d_m)/ir))
    p_err.append (t)
  noise = find_noise (p_err)
  f_pre_atrib = []
  j = 0
  for i in missing_list:
    f_pre_atrib.append(np.exp(noise[j] + yaf[i-1]))
    j = j+1
  return f_pre_atrib

def insert_dummyrows_withdensityvalues(cur, locid, noofrows,den_list):
  for i in xrange(0, noofrows):
    sql = "INSERT INTO retrospect_temp (MethodId,locid,density, tstrength_psi,tmodulus_psi,fstrength_psi,fmodulus_psi,hardness_inner,hardness_outer) VALUES ("+str(methid)+"," +str(locid) +","+str(den_list[i])+",null,null,null,null,null,null);"
    cur.execute(sql)
  cur.execute("commit")

def updatedb_withgenvals(cur, atr, genvals, sampleid_list):
  j = 0
  for item in sampleid_list: 
    # run the transformation
    sql = "UPDATE "+retro_cache.table_name +" SET " + str(atr) + " = " + str(genvals[j]) + " WHERE sampleid=" + str(item) + ";"
    cur.execute(sql)
    j = j + 1
  sql="commit"
  cur.execute(sql)


def generate_density_values(density_list, noofrows):
  n=0
  den_list = []
  minimum = int(min(density_list)*1000000)
  maximum = int(max(density_list)*1000000)
  for n in xrange(0,int(noofrows)):
    num=randrange(minimum,maximum)
    den_list.append(float(num)/float(1000000))
  return den_list


def find_markers_density(cur):
  sql = "SELECT tstrength_psi,sampleid FROM retrospect_temp" 
  cur.execute(sql)
  sampleid_list = []
  for row in cur.fetchall() :
    if row[0] == None:
      sampleid_list.append(row[1])
  return sampleid_list

def find_atrlist(cur,atr):
  atr_list = []
  sql = "select "+ str(atr)+ " from retrospect_temp where "+ str(atr)+ " is not NULL"
  cur.execute(sql)
  for row in cur.fetchall() :
    if row[0] != None:
      atr_list.append(row[0])
  return atr_list

#MAIN START
if ( len(sys.argv)  != 4 ):
        print 'Usage: retrospect <method id> <LocId> <NoofMockvalues> \n'
        print 'example: retrospect 1 10 23 \n Bye!!'
        sys.exit(2)

methid = sys.argv[1]
locat = sys.argv[2]
noofrows = sys.argv[3]
#code to generate random values for area

db = MySQLdb.connect(host="127.0.0.1", # your host, usually localhost
                     user="retrouser", # your username
                      passwd="retro123", # your password
                      db="retrodb") # name of the data base

# you must create a Cursor object. It will let
#  you execute all the query you need
cur = db.cursor()
 
sql = "truncate table retrospect_mock"
cur.execute(sql)
locations = []
if locat == '0':
  sql = "select locid from retrospect_locdb where flag = 1 and MethodId = " + str(methid)
  cur.execute (sql)
  for row in cur.fetchall():
     locations.append(row[0])
else:
  locations.append(locat)
print locations

for locat in locations:
  sql = "truncate table retrospect_temp"
  cur.execute(sql)
  cur.execute("commit")
  sql = "insert into retrospect_temp(MethodId,locid,density, tstrength_psi,tmodulus_psi,fstrength_psi,fmodulus_psi,hardness_inner,hardness_outer) SELECT MethodId,locid,density, tstrength_psi,tmodulus_psi,fstrength_psi,fmodulus_psi,hardness_inner,hardness_outer FROM retrospect WHERE locid =" + str(locat)
  cur.execute(sql)
  cur.execute("commit")
  atrib_list = ['tmodulus_psi','tstrength_psi','fstrength_psi','fmodulus_psi','hardness_inner','hardness_outer']
  retro_cache.table_name  = "retrospect_temp"


  for atr in atrib_list:
    sql = "SELECT density," + str(atr)+ ",sampleid FROM " + retro_cache.table_name +" where locid =" + str(locat)
    cur.execute(sql)
    atr_list = []
    density_list = []
    sampleid_list = []
    for row in cur.fetchall() :
      density_list.append(row[0])
      if row[1] == None:
        sampleid_list.append(row[2])
      atr_list.append(row[1])
    noofvals = 0
    noofvals = len(sampleid_list)
    if noofvals != 0:
      genvals = find_missingvalues(density_list, atr_list,sampleid_list)
      updatedb_withgenvals(cur, atr, genvals, sampleid_list)

  flag = 1
  for atr in atrib_list:
    if(flag == 1):
      den_list = generate_density_values( density_list, int(noofrows) )
      insert_dummyrows_withdensityvalues(cur, locat, int(noofrows),den_list)
      density_list.extend(den_list)
      samid_list = find_markers_density(cur)
      flag = flag + 1
    atr_list = find_atrlist(cur,atr)
    if len(samid_list) != 0:
      genvals = find_missingvalues(density_list, atr_list,samid_list)
      updatedb_withgenvals(cur, atr, genvals, samid_list)
  sql = "insert into retrospect_mock(MethodId,locid,density, tstrength_psi,tmodulus_psi,fstrength_psi,fmodulus_psi,hardness_inner,hardness_outer) SELECT MethodId,locid,density, tstrength_psi,tmodulus_psi,fstrength_psi,fmodulus_psi,hardness_inner,hardness_outer FROM retrospect_temp"
  cur.execute(sql)
  cur.execute("commit")

cur.close()

