# coding=utf-8
import numpy as np

import  math


class Joint():


    def __init__(self,index,type,length,theta = 0,offsets= None, frame = None):
        self.type = type
        self.frame = frame
        self.index = index
        self.linkln = length
        self.theta = theta
        self.offsets = offsets



class F_k():


    def get_proj_mat(self,joint1,joint2):
        c_theta = np.cos(joint1.theta)
        f1 = joint1.frame
        f2 = joint2.frame
        final_proj_mat = None
        for i in range(0,3):
            pro_mt = [[], [], []]
            if (f2[i][0] == 1):
                a = 1
            else:
                a = -1
            if(f1[i][1] == 'x'):
                if (f2[i][1] == 'x'):
                    pro_mt = np.matrix([[a],[0],[0]])
                elif (f2[i][1] == 'y'):
                    pro_mt = np.matrix([[0],[a],[0]])
                elif (f2[i][1] == 'z'):
                    pro_mt = np.matrix([[0],[0],[a]])
            elif(f1[i][1]=='y'):
                if (f2[i][1] == 'x'):
                    pro_mt = np.matrix([[a],[0],[0]])
                elif (f2[i][1] == 'y'):
                    pro_mt = np.matrix([[0],[a],[0]])
                elif (f2[i][1] == 'z'):
                    pro_mt = np.matrix([[0],[0],[a]])
            elif (f1[i][1] == 'z'):
                if (f2[i][1] == 'x'):
                    pro_mt = np.matrix([[a],[0],[0]])
                elif (f2[i][1] == 'y'):
                    pro_mt = np.matrix([[0],[a],[0]])
                elif (f2[i][1] == 'z'):
                    pro_mt = np.matrix([[0],[0],[a]])
            final_proj_mat = np.concatenate((final_proj_mat, pro_mt), 1)
        return final_proj_mat


    def cal_disp_mat(self, disp, offsets, type):
        if type =='revolute':
            cos = np.cos((disp*np.pi/180))
            sin = np.sin((disp*np.pi/180))
            return np.matrix([[offsets[0] * cos], [offsets[1] * sin], [offsets[2]]])
        elif type=='prismatic':
            return np.matrix([[offsets[0]], [offsets[1]], [offsets[2] + disp]])

    def get_Z_RM(self,theta,type):
        if type == "revolute":
            cos = np.cos((theta * np.pi / 180))
            sin = np.sin((theta * np.pi / 180))
            return np.matrix([[cos,sin,0],[-sin,cos,0],[0,0,1]])
        elif type=="prismatic":
            return np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def final_rtmat(self,proj_mat,rot_mat):
        return np.multiply(rot_mat,proj_mat)

    def get_H_mat(self,fnl_rot_mat,disp_mat):
        h_mat= np.concatenate((fnl_rot_mat,disp_mat),1)
        return np.concatenate((h_mat,np.matrix([[0,0,0,1]])),0)

class DH_():

    theta = []
    alpha = []
    r = []
    d = []
    HM_ls = []
    HM = 1
    def __init__(self,alpha,r,d):# in the table
        self.alpha = alpha
        self.r = r
        self.d=d

    def __calc_H_Matrix__(self):
        if(len(self.theta) == len(self.alpha) == len(self.r) == len(self.d)):
            for i in range(0,len(self.theta)-1) :
                theta = (self.theta[i])*np.pi/180.0;alpha = (self.alpha[i])*np.pi/180.0
                r=self.r[i];d=self.d[i]
                st = np.sin(theta);ct=np.cos(theta);sa=np.sin(alpha);ca=np.cos(alpha)
                H = [[ct,-st*ca,st*sa,r*ct],[st,ct*ca,-ct*sa,r*st],[0,sa,ca,d],[0,0,0,1]]
                self.HM_ls.append(H)
                self.HM=np.multiply(self.HM,H)
        else:
             IndexError.message("length of parameters is not same")
    def show_HM(self):
        return self.HM
    def calc_effector_Pos(self,theta):
        self.theta = theta
        self.__calc_H_Matrix__()
        return self.HM




class Arm():

    joints = []
    joints.append(Joint(0, 0, 0, 0, 0))
    implementation = {}
    def __init__(self,approach):
        if approach == "forward kinematics":
            proj_mats = []
            disp_mats = []
            f = F_k()
            self.implementation["approach"] = "F_K"
            self.implementation["object"] = F_k()
            self.implementation["proj_mats"] = proj_mats

    def __Add_joint__(self, type, length, theta=0, offsets=None, frame=None):
        self.joints.append(Joint(len(self.joints), type, length, theta, offsets, frame))
        if self.implementation["approach"] == "F_K":
            if len(self.joints) > 1:
                self.implementation["proj_mats"].append(
                    self.implementation["object"].get_proj_mat(self.joints[len(self.joints) - 2],
                                                               self.joints[len(self.joints) - 2]))

    def calc_change(self,disp_vec):
        f = F_k()
        if(self.implementation["approach"]=="F_K"):     #format disp_vec = [j1_disp/theta,-----]
            c = 0
            H_mt = []
            for i in disp_vec:
                Z_Rot_mat = f.get_Z_RM(i,self.joints[c].type)
                Prj_mat = self.implementation["proj_mats"][c]
                Rot_mat = np.multiply(Z_Rot_mat,Prj_mat)
                dis_mat = f.cal_disp_mat(i,self.joints[c].offsets,self.joints[c].type)
                H_mt.append(self.implementation["object"].get_H_mat(Rot_mat,dis_mat))
            H = np.identity(3)
            for i in range(0,len(H_mt)):
                H = np.multiply(H,H_mt[i])
            return H




