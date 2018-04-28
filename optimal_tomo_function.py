# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 14:15:43 2016

@author: jmcui
"""
#from TomographyQubitCoolingValidity import *
import scipy.linalg as linalg
import numpy as np
from scipy.optimize import minimize, fmin_cg, fmin


def estimate_function_old(pmt, prob, ta1_initial):  # ¸ø¶¨¸ÅÂÊ£¬¸ø³ö²ÎÊýº¯Êý
    premeter = np.concatenate(([ta1_initial], pmt))  # È·¶¨t1µÄÖµ£¬ÒòÎª¾Å¸ö²ÎÊýÖ»ÓÐ°Ë¸ö±äÁ¿
    return estimate_old(premeter, prob)

def estimate_function(pmt, prob, ta1_initial):  # ¸ø¶¨¸ÅÂÊ£¬¸ø³ö²ÎÊýº¯Êý
    premeter = np.concatenate(([ta1_initial], pmt))  # È·¶¨t1µÄÖµ£¬ÒòÎª¾Å¸ö²ÎÊýÖ»ÓÐ°Ë¸ö±äÁ¿
    return estimate(premeter, prob)
    
def estimate_old(t, ni):  # t为待估计的参数，ni为样本次数，而此时总数为1，所以ni为频率
    '''# 这个函数对应8个正常tomo的测量结果
    '''
    r_t = np.linalg.norm(t)**2  # 求TT的迹
    a1 = t[0]**2 + t[3]**2 + t[4]**2 + t[7]**2 + t[8]**2
    a2 = t[2]**2
    a3 = t[0]**2 / 2.0 + t[1]**2 / 2.0 + t[1] * t[3] + t[3]**2 / 2.0 + t[4]**2 / 2.0 + \
        t[5]**2 / 2.0 + t[5] * t[7] + t[6]**2 / 2.0 + t[6] * t[8] + t[7]**2 / 2.0 + t[8]**2 / 2.0
    a4 = t[0]**2 / 2.0 + t[1]**2 / 2.0 + t[1] * t[4] + t[3]**2 / 2.0 + t[4]**2 / 2.0 + \
        t[5]**2 / 2.0 + t[5] * t[8] + t[6]**2 / 2.0 - t[6] * t[7] + t[7]**2 / 2.0 + t[8]**2 / 2.0
    a5 = t[0]**2 / 2.0 + t[2]**2 / 2.0 + t[2] * t[7] + t[3]**2 / 2.0 + t[4]**2 / 2.0 + t[7]**2 / 2.0 + t[8]**2 / 2.0
    a6 = t[0]**2 / 2.0 + t[2]**2 / 2.0 + t[2] * t[8] + t[3]**2 / 2.0 + t[4]**2 / 2.0 + t[7]**2 / 2.0 + t[8]**2 / 2.0
    a7 = t[1]**2 / 2.0 + t[2]**2 / 2.0 + t[2] * t[5] + t[5]**2 / 2.0 + t[6]**2 / 2.0
    a8 = t[1]**2 / 2.0 + t[2]**2 / 2.0 + t[2] * t[6] + t[5]**2 / 2.0 + t[6]**2 / 2.0
    b = np.array([a1, a2, a3, a4, a5, a6, a7, a8])
    # print np.sum(ni*np.log(b/r_t))
    # 负号是因为我们利用的是求最小值函数fmin，得到最大值结果
    return -1.0 * np.sum(ni * np.log(b / r_t) + (1.0 - ni) * np.log(1.0 - b / r_t))


def estimate(t, ni):  # t为待估计的参数，ni为样本次数，而此时总数为1，所以ni为频率
    '''这个函数对应12个无偏正交基的测量结果
    '''
    r_t = np.linalg.norm(t)**2  # 求TT的迹
    c = 1 / 3.0
    c1 = t[0]**2 + t[1]**2 + t[2]**2 + t[3]**2 + t[4]**2 + t[5]**2 + t[6]**2 + t[7]**2 + t[8]**2
    c2 = t[1] * t[3] + t[2] * t[5] + t[2] * t[7] + t[5] * t[7] + t[6] * t[8]
    a1 = c1 + 2 * c2
    a2 = c1 - c2 + np.sqrt(3.0) * (t[1] * t[4] + t[2] * t[6] - t[2] * t[8] + t[5] * t[8] - t[6] * t[7])
    a3 = c1 - c2 - np.sqrt(3.0) * (t[1] * t[4] + t[2] * t[6] - t[2] * t[8] + t[5] * t[8] - t[6] * t[7])
    a4 = 3 * (t[0]**2 + t[3]**2 + t[4]**2 + t[7]**2 + t[8]**2)
    a5 = 3 * (t[1]**2 + t[5]**2 + t[6]**2)
    a6 = 3 * t[2]**2
    a7 = c1 + (2 * t[1] * t[3] - t[2] * t[5] - t[2] * t[7] + 2 * t[5] * t[7] +
               2 * t[6] * t[8]) + np.sqrt(3.0) * (-t[2] * t[6] - t[2] * t[8])
    a8 = c1 + (-t[1] * t[3] + 2 * t[2] * t[5] - t[2] * t[7] - t[5] * t[7] - t[6] * t[8]) + \
        np.sqrt(3.0) * (t[1] * t[4] + t[2] * t[8] + t[5] * t[8] - t[6] * t[7])
    a9 = c1 + (-t[1] * t[3] - t[2] * t[5] + 2 * t[2] * t[7] - t[5] * t[7] - t[6] * t[8]) + \
        np.sqrt(3.0) * (-t[1] * t[4] + t[2] * t[6] - t[5] * t[8] + t[6] * t[7])
    a10 = c1 + 2 * t[1] * t[3] - t[2] * t[5] - t[2] * t[7] + 2 * t[5] * \
        t[7] + 2 * t[6] * t[8] + np.sqrt(3.0) * (t[2] * t[6] + t[2] * t[8])
    a11 = c1 + (-t[1] * t[3] - t[2] * t[5] + 2 * t[2] * t[7] - t[5] * t[7] - t[6] * t[8]) - \
        np.sqrt(3.0) * (-t[1] * t[4] + t[2] * t[6] - t[5] * t[8] + t[6] * t[7])
    a12 = c1 + (-t[1] * t[3] + 2 * t[2] * t[5] - t[2] * t[7] - t[5] * t[7] - t[6] * t[8]) - \
        np.sqrt(3.0) * (t[1] * t[4] + t[2] * t[8] + t[5] * t[8] - t[6] * t[7])

    b = np.array([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12]) * c
    # print np.sum(ni*np.log(b/r_t))
    return -1.0 * np.sum(ni * np.log(b / r_t) + (1.0 - ni) * np.log(1.0 - b / r_t))  # 负号是因为我们利用的是求最小值函数fmin，得到最大值结果


tes = np.array([0.19333333, 0.49, 0.26833333, 0.425, 0.47666667, 0.605, 0.73, 0.13])
def optimal_rho_old(probi):  # Íâ²¿µ÷ÓÃÕâ¸öº¯Êý
    ta1_initial = 0.0
    t0 = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])  # initial values of parameters
    ta = fmin(lambda x: estimate_function_old(x, probi, ta1_initial), t0,
              disp=False)  # Çó²ÎÊýº¯ÊýµÄ×î´óÖµÊ±µÄ²ÎÊýÖµ£¬È»ºóÖØ¹¹¾ØÕó
    ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8, ta9 = np.concatenate(([ta1_initial], ta))
    T_estimate = np.array([[ta1, ta4 + 1.0j * ta5, ta8 + 1.0j * ta9], [0, ta2, ta6 + 1.0j * ta7], [0, 0, ta3]])
    M = T_estimate.dot(T_estimate.conjugate().T)
    M = M / np.trace(M)
    return M

def optimal_rho(probi):  # Íâ²¿µ÷ÓÃÕâ¸öº¯Êý
    ta1_initial = 0.0
    t0 = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])  # initial values of parameters
    ta = fmin(lambda x: estimate_function(x, probi, ta1_initial), t0,
              disp=False)  # Çó²ÎÊýº¯ÊýµÄ×î´óÖµÊ±µÄ²ÎÊýÖµ£¬È»ºóÖØ¹¹¾ØÕó
    ta1, ta2, ta3, ta4, ta5, ta6, ta7, ta8, ta9 = np.concatenate(([ta1_initial], ta))
    T_estimate = np.array([[ta1, ta4 + 1.0j * ta5, ta8 + 1.0j * ta9], [0, ta2, ta6 + 1.0j * ta7], [0, 0, ta3]])
    M = T_estimate.dot(T_estimate.conjugate().T)
    M = M / np.trace(M)
    return M
