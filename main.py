#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from osgeo import ogr
from osgeo import osr
import math as m
import pyproj


def wgs_to_jtsk(point):
    """
    Transformation of point geometry from EPSG:4326 to EPSG:5514
    :param point: point geometry in EPSG:4326
    :return: X and Y coordinate in EPSG:5514
    """
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)
    target = osr.SpatialReference()
    target.ImportFromEPSG(5514)
    transform = osr.CoordinateTransformation(source, target)
    point.Transform(transform)
    point.ExportToWkt()
    return point.GetX(), point.GetY()


def create_point(lat, long):
    point_wkt = "POINT (%s %s)" % (lat, long)
    return ogr.CreateGeometryFromWkt(point_wkt)


def traffic_sign_panos(dataframe):
    """
    Creating list of panoramic images with the same traffic sign.
    :return: List of panoramic images with the same
    """
    pass


def traffic_sign_azimuth(centroid):
    """
    Calculating azimuth of traffic sign based on its position in panoramic image.
    :param centroid: list of centroids of bounding boxes
    :return:
    """

    pass


def traffic_sign_table(pano_table, traffic_sign_pairs):
    M_id = traffic_sign_pairs['M_id']
    N_id = traffic_sign_pairs['N_id']
    pano_ids = sorted(list(set(M_id) | set(N_id)))
    pano_subset = pano_table[pano_table['panorama_file_name'].isin(pano_ids)]
    return pano_subset


def traffic_sign_location(row):
    M = [row['M_lat'], row['M_long']]
    N = [row['N_lat'], row['N_long']]
    dfi_MN = N[0] - M[0]
    dlambda_MN = N[1] - M[1]

    A_Mj = (row['M_heading'] + 180) / 180 * m.pi
    A_Nj = (row['N_heading'] + 180) / 180 * m.pi
    # pix_M = 686 / 8000 * 360
    # pix_N = 1200 / 8000 * 360
    pix_M = row['M_X'] / 8000 * 360
    pix_N = row['N_X'] / 8000 * 360
    A_Mz = pix_M * m.pi / 180
    A_Nz = pix_N * m.pi / 180
    A_Mc = A_Mj + A_Mz
    A_Nc = A_Nj + A_Nz

    #A_M_roll = row['M_roll']
    #A_M_Y = (row['M_Y'] / 2000 * 90)
    #A_M_yc = (A_M_roll + A_M_Y) / 180 * m.pi

    d = m.sin(A_Mc - A_Nc)
    S_AP = 1 / d * (m.cos(A_Nc) * dlambda_MN - m.sin(A_Nc) * dfi_MN)
    #S_znacka = S_AP/m.cos(A_M_yc)
    return M[0] + S_AP * m.cos(A_Mc), M[1] + S_AP * m.sin(A_Mc)


def traffic_sign_loc(row):
    # load points
    m_x = row['M_long']
    m_y = row['M_lat']
    n_x = row['N_long']
    n_y = row['N_lat']

    # load headings
    m_h = m.radians(row['M_heading'])
    n_h = m.radians(row['N_heading'])

    if m_h < m.pi:
        m_head = m_h + m.pi
    else:
        m_head = m_h - m.pi

    if n_h < m.pi:
        n_head = n_h + m.pi
    else:
        n_head = n_h - m.pi

    print("m_head: ", m.degrees(m_head))
    print("n_head: ", m.degrees(n_head))

    # calculate delta
    dx = n_x - m_x
    dy = n_y - m_y
    d = m.atan(dx / dy)
    if dx > 0 < dy:
        delta = d
    elif dx > 0 > dy:
        delta = m.pi - d
    elif dx < 0 > dy:
        delta = m.pi + d
    elif dx < 0 < dy:
        delta = 2*m.pi - d

    print("dx: ", dx)
    print("dy: ", dy)
    print("delta: ", m.degrees(delta))

    # calculate pix angles
    m_pix = m.radians(row['M_X'] / 8000 * 360)
    n_pix = m.radians(row['N_X'] / 8000 * 360)

    print("m_pix: ", m.degrees(m_pix))
    print("n_pix: ", m.degrees(n_pix))

    # calculate sigma deviations

    if delta - m_head > 0:
        m_sig = delta - m_head
    else:
        m_sig = m.pi + (delta - m_head)

    if delta - n_head > 0:
        n_sig = delta - n_head
    else:
        n_sig = m.pi + (delta - n_head)

    print("m_sig: ", m.degrees(m_sig))
    print("n_sig: ", m.degrees(n_sig))

    # calculate alpha and beta

    if (m_pix - m_sig) + (n_pix - n_sig) < m.pi:
        alpha = m_pix - m_sig
        beta = n_pix - n_sig
    else:
        alpha = m.pi - (m_pix - m_sig)


    if n_pix - n_sig < m.pi:
        beta = n_pix - n_sig
    else:
        beta = m.pi - (n_pix - n_sig)

    print("alpha: ", m.degrees(alpha))
    print("beta: ", m.degrees(beta))

    # calculate coordinates of P
    a_cot = 1 / m.tan(alpha)
    b_cot = 1 / m.tan(beta)
    q = a_cot + b_cot
    p_x = ((n_y - m_y) + n_x * b_cot + m_x * a_cot) / q
    p_y = ((m_x - n_x) + n_y * b_cot + m_y * a_cot) / q

    return p_x, p_y


def traffic_sign_loc2(row):
    # load points
    m_x = row['M_long']
    m_y = row['M_lat']
    n_x = row['N_long']
    n_y = row['N_lat']

    print("M: ", m_y, m_x)
    print("N: ", n_y, n_x)

    # load headings
    m_head = (row['M_heading'] + 180) % 360
    n_head = (row['N_heading'] + 180) % 360

    print("M_heading: ", m_head)
    print("N_heading: ", n_head)

    # azimuth
    geodesic = pyproj.Geod(ellps='WGS84')
    m_azimuth, n_azimuth, distance = geodesic.inv(m_x, m_y, n_x, n_y)
    m_azimuth = m_azimuth % 360
    n_azimuth = n_azimuth % 360

    print("fwd_azimuth: ", m_azimuth)
    print("back_azimuth: ", n_azimuth)

    # pix angles
    m_pix = row['M_X'] / 8000 * 360
    n_pix = row['N_X'] / 8000 * 360

    print("m_pix: ", m_pix)
    print("n_pix: ", n_pix)

    # dir angles
    m_dir = (m_head + m_pix) % 360
    n_dir = (n_head + n_pix) % 360

    print("m_dir: ", m_dir)
    print("n_dir: ", n_dir)

    # alpha, beta
    alpha = m.radians(abs(m_dir - m_azimuth) % 360)
    beta = m.radians(abs(n_dir - n_azimuth) % 360)

    print("alpha: ", m.degrees(alpha))
    print("beta: ", m.degrees(beta))

    # calculate coordinates of P
    a_cot = 1 / m.tan(alpha)
    b_cot = 1 / m.tan(beta)
    q = a_cot + b_cot
    p_x = ((n_y - m_y) + n_x * b_cot + m_x * a_cot) / q
    p_y = ((m_x - n_x) + n_y * b_cot + m_y * a_cot) / q

    return p_y, p_x


def traffic_sign_loc3(row):
    # load points
    m_x = row['M_long']
    m_y = row['M_lat']
    n_x = row['N_long']
    n_y = row['N_lat']

    print("M: ", m_y, m_x)
    print("N: ", n_y, n_x)

    # load headings
    m_head = (row['M_heading'] + 180) % 360
    n_head = (row['N_heading'] + 180) % 360

    print("M_heading: ", m_head)
    print("N_heading: ", n_head)

    # pix angles
    m_pix = row['M_X'] / 8000 * 360
    n_pix = row['N_X'] / 8000 * 360

    print("m_pix: ", m_pix)
    print("n_pix: ", n_pix)

    # dir angles
    m_dir = m.radians((m_head + m_pix) % 360)
    n_dir = m.radians((n_head + n_pix) % 360)

    print("m_dir: ", m_dir)
    print("n_dir: ", n_dir)

    q1 = m.tan(m_dir) - m.tan(n_dir)
    cot_a = 1/m.tan(m_dir)
    cot_b = 1/m.tan(n_dir)
    q2 = cot_a - cot_b
    p_y = (n_x - m_x + m_y * m.tan(m_dir) - n_y * m.tan(n_dir)) / q1
    p_x = (n_y - m_y + m_x * cot_a - n_x * cot_b) / q2

    return p_y, p_x


# Calculates Rotation Matrix given euler angles.
def rotation_matrix(omega, phi, kappa):
    R_x = np.array([[1, 0, 0],
                    [0, m.cos(omega), -m.sin(omega)],
                    [0, m.sin(omega), m.cos(omega)]
                    ])

    R_y = np.array([[m.cos(phi), 0, m.sin(phi)],
                    [0, 1, 0],
                    [-m.sin(phi), 0, m.cos(phi)]
                    ])

    R_z = np.array([[m.cos(kappa), -m.sin(kappa), 0],
                    [m.sin(kappa), m.cos(kappa), 0],
                    [0, 0, 1]
                    ])

    R = np.dot(R_z, np.dot(R_y, R_x))
    return R


def polar2world(lat, lon):
    x = m.sin(lat) * m.cos(lon)
    y = m.sin(lat) * m.sin(lon)
    z = m.cos(lat)
    return np.array([x, y, z])


def world2polar(arr):
    lat = m.acos(arr[2])
    lon = m.asin(arr[1]/m.sin(lat))
    return lat, lon


def traffic_sign_loc4(row):
    # GPS locations of panoramic images
    m_lat = row['M_lat']
    m_lon = row['M_long']
    m_el = row['M_el']
    n_lon = row['N_long']
    n_lat = row['N_lat']
    n_el = row['N_el']

    # image coordinates of object
    m_phi = m.radians(row['M_X'] / 4000 * 180)
    m_lam = m.radians(row['M_Y'] / 8000 * 360)
    n_phi = m.radians(row['N_X'] / 4000 * 180)
    n_lam = m.radians(row['N_Y'] / 8000 * 360)

    # rotation angles of panoramic images
    m_roll = m.radians(row['M_roll'])
    m_pitch = m.radians(row['M_pitch'])
    m_yaw = m.radians(row['M_heading'])
    n_roll = m.radians(row['N_roll'])
    n_pitch = m.radians(row['N_pitch'])
    n_yaw = m.radians(row['N_heading'])

    # phi,lam -> X, Y, Z
    m_arr = polar2world(m_phi, m_lam)
    n_arr = polar2world(n_phi, n_lam)

    # rotation matrices
    m_rotation = rotation_matrix(m_roll, m_pitch, m_yaw)
    n_rotation = rotation_matrix(n_roll, n_pitch, n_yaw)

    # rectification of XYZ
    m_arr_rec = np.invert(m_rotation) * np.transpose(m_arr)
    n_arr_rec = np.invert(n_rotation) * np.transpose(n_arr)

    # XYZ -> phi,lam
    m_lat2, m_lon2 = world2polar(np.transpose(m_arr_rec))
    n_lat2, n_lon2 = world2polar(np.transpose(n_arr_rec))

    print("new M direction: ", m.degrees(m_lat2), m.degrees(m_lon2))
    print("new N direction: ", m.degrees(n_lat2), m.degrees(n_lon2))


def traffic_signs_shapefile(path, table):
    driver = ogr.GetDriverByName('ESRI Shapefile')  # set driver and vector format
    if os.path.exists(path):
        driver.DeleteDataSource(path)
    datasource = driver.CreateDataSource(path)  # create datasource
    layer = datasource.CreateLayer('points', geom_type=ogr.wkbPoint)  # create new layer linestring
    # Layer is still empty (but has two attribute columns), now add features
    for line in table:
        featureDefn = layer.GetLayerDefn()  # define new feature on Layer
        feature = ogr.Feature(featureDefn)  # create new feature
        point = ogr.Geometry(ogr.wkbPoint)  # generate point geometry
        point.AddPoint(line[1], line[2])  # add vertices to linestring
        feature.AddGeometry(point)  # set geometry of feature using linestring
        layer.CreateFeature(feature)  # create feature and write to layer
        feature = None  # clean up
    datasource = None


def table_to_shp(path, csv):
    csv_gdf = gpd.GeoDataFrame(csv, geometry=gpd.points_from_xy(csv.lat, csv.long), crs='epsg:4326')
    csv_gdf.to_file(path)


if __name__ == '__main__':
    # table_pano_path = r'E:\Gekon\znacky\data\Praha21Bechexp_panorama.csv'
    # table_pano = pd.read_csv(table_pano_path, delimiter=';', usecols=[1, 2, 3, 7])
    # print(table_pano)
    table_znacky_path = r'C:\Users\micha\Desktop\Gekon\znacky\data\praha_1\znacky_new2.csv'
    table_znacky = pd.read_csv(table_znacky_path, delimiter=';')
    table_znacky.apply(pd.to_numeric, errors='coerce').fillna(table_znacky)

    pano_table_path = r'C:\Users\micha\Desktop\Gekon\znacky\data\praha_1\PrahaCentrum1exp_panorama.csv'
    pano_table = pd.read_csv(table_znacky_path, delimiter=';')

    traffic_sign_pairs_path = r'C:\Users\micha\Desktop\Gekon\znacky\data\praha_1\znacky_list.csv'
    traffic_sign_pairs = pd.read_csv(traffic_sign_pairs_path, delimiter=';')



    for i in range(11):
        lat, long = traffic_sign_loc3(table_znacky.iloc[i])
        table_znacky.loc[i, 'lat'] = float(lat)
        table_znacky.loc[i, 'long'] = float(long)
        print(lat, long)

    table_znacky.to_csv(r'C:\Users\micha\Desktop\Gekon\znacky\data\praha_1\znacky_new3.csv', sep=';')
    table_to_shp(r'C:\Users\micha\Desktop\Gekon\znacky\data\praha_1\table_znacky.shp', table_znacky)
    # long, lat = traffic_sign_location(table_pano,
    # P_point = create_point(long, lat)
    # X, Y = wgs_to_jtsk(P_point)
