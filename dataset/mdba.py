# -*- coding:UTF-8 -*-
# Multiobjective Discrete Bat Algorithm 多目标离散蝙蝠算法
"""
Input:A dynamic network DN = {G1,G2,...,GT}
Output:Network partition at each time P = {P1,P2,...,PT}
"""
import networkx as nx
import util.graph_helper as gh
import right_partition2 as rp
import util.modularity as mq
from util.nmi_test import NMI
from util.nmi_test2 import NMI2
from util.common_nodes import cal_com_node
from dba import dba
from util.pareto import pareto
from util.error_rate import com_error_rate
import random
import time
# 初始化算法的基础参数：种群大小N（即蝙蝠数量：每个蝙蝠代表一个划分结果），
# 时间步T的数量，代数的最大值genmax（迭代次数）


def calc_par(num):
    T = 11  # 总时间步 一共15个时间步
    bat_N = 100  # 蝙蝠个数
    # 使用CNM算法获得图G1的分区结果P1 = {p11,p12,...,P1m}(1下标表示时间片，m上标表示社区)
    # 此处见CNM.py
    # data/syntetic2/1.edgelist
    # data2/synfix/z_3/synfix_3.t
    # data/2/switch.t01.comm
    # data2/synvar/z_3/synvar_3.t01.comm1
    # data/Cell/real.t0
    file_qian = "data2/synfix/z_5/synfix_5.t"
    file_she = ".comm"
    file_bian = ".edges"
    G_pre = gh.load_graph(file_qian + "01" + file_bian)
    G_pre_nn = G_pre.number_of_nodes()  # 前一个时间片图的节点个数
    components_pre = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64], [65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96], [128, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127]]
    # print components_pre
    # 1到T时间步
    # [pop[bat]]
    for t in range(2, T):  # (1,T)
        # 开始时间
        start = time.clock()
        # 加载数据
        """注意引用文件名"""
        if t < 10:
            t = "0" + str(t)
        G = gh.load_graph(file_qian + str(t) + file_bian)
        G_nn = G.number_of_nodes()  # 图G的节点个数
        com_nn = cal_com_node(G_pre, G)
        # 蝙蝠种群集合，第i个蝙蝠,一个蝙蝠代表一个解，一个解有d维，即节点的个数，列表位置为节点，结果为节点随机邻居
        pop = []  # 邻接表示结合
        pop_partition = []  # 蝙蝠种群集合，分区集合
        pop_value = []  # 分区多目标结果集合
        rep = []  # 帕累托最优解邻接表示集合
        rep_partition = []  # 分区结合
        # 初始化蝙蝠集合，即随机N个情况
        for i in range(bat_N):  # (0, bat_N)
            # G2为解码后的图
            G2 = nx.Graph()
            # 基于轨迹的邻接表示
            pop_bat = []  # 初始化每个蝙蝠
            for node in G.nodes():
                # print j, [n for n in G.neighbors(j)]
                random_nb = random.sample([n for n in G.neighbors(node)], 1)[0]
                pop_bat.append([node,random_nb])
                G2.add_edge(node, random_nb)
            # print pop_bat
            pop.append(pop_bat)  # 蝙蝠种群集合pop
            # 不相交集算法
            components = [list(c) for c in list(nx.connected_components(G2))]
            # print "社区为："+str(components)
            Q = mq.cal_Q(components, G)
            # 每次与前一个进行比较
            # a = NMI2(components_pre, components, G_pre_nn,G_nn,com_nn)
            a = NMI(components_pre, components, G_nn)
            nmi = a.nmi()
            # print Q
            # print nmi
            pop_partition.append(components)  # 蝙蝠分区情况集合
            pop_value.append([i, Q, nmi])  # 蝙蝠目标值集合
        # print len(pop), len(result_mul)
        # print pop
        # print pop_partition
        # print pop_value
        # 帕累托最优解集合
        rep_value = pareto(pop_value)
        print "rep_value:",rep_value
        for i in range(len(rep_value)):
            j = rep_value[i][0]
            rep.append(pop[j])
            rep_partition.append(pop_partition[j])
        # print rep
        # print rep_partition

        # DBA算法
        # 选择一个蝙蝠
        # 解码蝙蝠获得时间步t时刻的分区
        best = dba(G, pop, pop_value, pop_partition, rep, rep_value, rep_partition, components_pre,G_nn,G_pre_nn,com_nn)  # G_nn,G_pre_nn,com_nn
        print t, best[0], best[1]

        # 计算错误率
        components_right = rp.calc(file_qian + str(t) + file_she)
        print "正确划分：", components_right
        error = com_error_rate(components_right, best[1])
        print "error:", error
        # 计算与标准数据NMI值
        b = NMI(components_right, best[1], G_nn)
        nmi_gt = b.nmi()
        print "与标准结果NMI：", nmi_gt

        # 结束时间
        elapsed = (time.clock() - start)
        print "Time used in",t,":",elapsed
        print "--------------------------------------------------------------------"
    return 0


if __name__ == "__main__":
    calc_par(1)


