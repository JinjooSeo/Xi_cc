#include <cmath>
#include "AliDummyHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisTaskSEXiccTopKpipi.h"
#include "TStopwatch.h"
// using namespace std;

//SETTINGS
//************************************

Bool_t runLocal=kTRUE;                                  // flag to run locally on AliAOD.root + AliAOD.VertexingHF.root
//TString pathToLocalAODfiles="../Omega_ccc/input/bkg/";//"./analysis/background"; // ../input_files/bkg_Kr path to find AOD files when running locally
TString pathToLocalAODfiles="../data/210115";
Bool_t runGridTest=kFALSE;                                // flag to run a grid test: kTRUE (+runLocal=kFALSE). To run job on GRID: runGridTest=kFALSE, runLocal=kFALSE
TString runMode="full";                                  // sets the run grid mode: "full", "terminate"

TString aliPhysVersion="vAN-20200721_ROOT6-1";

Bool_t isRunOnMC=kTRUE;                                 // set to kTRUE to run on Mone Carlo and uncomment/comment accordingly the following lines about paths on Alien
//paths on Alien
// Monte Carlo
TString gridDataDir="/alice/cern.ch/user/a/afestant/testHIJINGPYTHIAgen_omegaccc_bkg/output_Kr";//"/alice/sim/2016/LHC16i2a/";
//TString gridDataDir="/alice/cern.ch/user/s/strogolo/Background_Omega_ccc";
TString gridDataPattern="";//"/AOD198";
// Data
//TString gridDataDir="/alice/data/2015/LHC15o/";
//TString gridDataPattern="/pass1_pidfix/AOD194";


// Alien output directory
TString gridWorkingDir="testHIJINGPYTHIAgen_omegac_bkg/config_0/background_analysis_full_Kr_4th";
TString gridOutputDir="output";

//run numbers
const Int_t nruns = 130;
//Int_t runlist[nruns] = {2};
//Int_t runlist[nruns] = {245146};
//Int_t runlist[nruns] = {1, 2, 3, 4, 5, 6, 7, 8, 9};//nruns = 9
/*Int_t runlist[nruns] = {
  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
  26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
  42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
  58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
  74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
  90, 91, 92, 93, 94, 95, 96, 97, 98, 99
};*/ //nruns = 90
Int_t runlist[nruns] = {
  /*100, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, //11
  101, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 102, 1020, //13
  1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 103, 1030, 1031, 1032,
  1033, 1034, 1035, 1036, 1037, 1038, 1039, 104, 1040, 1041, 1042, 1043, 1044,
  1045, 1046, 1047, 1048, 1049, 105, 1050, 1051, 1052, 1053, 1054, 1055, 1056,
  1057, 1058, 1059, 106, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068,
  1069, 107, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 108,
  1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 109, 1090, 1091};*/
  1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 110, 1100, 1101, 1102, 1103,
  1104, 1105, 1106, 1107, 1108, 1109, 111, 1110, 1111, 1112, 1113, 1114, 1115,
  1116, 1117, 1118, 1119, 112, 1120, 1121, 1122, 1123, 1124, 1125, 1127, 1128,
  1129, 113, 1130, 1131, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 114, 1140,
  1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 115, 1150, 1151, 1152,
  1153, 1154, 1155, 1156, 1157, 1158, 1159, 116, 1160, 1161, 1162, 1163, 1164,
  1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177,
  1178, 1179, 118, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189,
  119, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199, 120, 1200,
  1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 121, 1210, 1211, 1212};
  /*1213, 1214, 1215, 1216, 1217, 1218, 1219, 122, 1220, 1221, 1222, 1223, 1224,
  1225, 1226, 1227, 1228, 1229, 123, 1230, 1231, 1232, 1233, 1234, 1235, 1236,
  1237, 1238, 1239, 124, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248,
  1249, 125, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 126,
  1260, 1261, 1262, 1263, 1264, 1265, 1266, 1267, 1268, 1269, 127, 1270, 1271,
  1272, 1273, 1274, 1275, 1276, 1277, 1279, 128, 1280, 1281, 1282, 1283, 1284,
  1285, 1286, 1287, 1288, 1289, 129, 1290, 1291, 1292, 1293, 1294, 1295, 1296,
  1297, 1298, 1299, 130, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308,
  131, 1310, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319, 132, 1320,
  1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 133, 1330, 1331, 1332,
  1333, 1334, 1335, 1336, 1337, 1338, 1339, 134, 1340, 1341, 1342, 1343, 1344,
  1345, 1346, 1347, 1348, 1349, 135, 1350, 1351, 1352, 1353, 1354, 1355, 1356,
  1357, 1358, 1359, 136, 1360, 1361, 1362, 1363, 1364, 1365, 1366, 1367, 1368,
  1369, 137, 1370, 1371, 1372, 1373, 1374, 1375, 1376, 1377, 1378, 1379, 138,
  1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 139, 1390, 1391,
  1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 140, 1401, 1402, 1403, 1404,
  1405, 1406, 1407, 1408, 1409, 141, 1410, 1411, 1412, 1413, 1414, 1415, 1416,
  1417, 1418, 1419, 142, 1420, 1421, 1422, 1423, 1424, 1425, 1426, 1427, 1428,
  1429, 143, 1430, 1431, 1432, 1433, 1434, 1435, 1436, 1437, 1438, 1439, 144,
  1440, 1441, 1442, 1443, 1444, 1445, 1446, 1447, 1448, 1449, 145, 1450, 1451,
  1452, 1453, 1454, 1455, 1456, 1457, 1458, 1459, 146, 1460, 1461, 1462, 1463,
  1464, 1465, 1466, 1467, 1468, 1469, 147, 1470, 1471, 1472, 1473, 1474, 1475,
  1476, 1477, 1478, 1479, 148, 1480, 1481, 1482, 1483, 1484, 1485, 1486, 1487,
  1488, 1489, 149, 1490, 1491, 1492, 1493, 1494, 1495, 1496, 1497, 1498, 150,
  1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508, 1509, 151, 1510, 1511,
  1512, 1513, 1514, 1515, 1516, 1517, 1518, 152, 1520, 1521, 1522, 1523, 1524,
  1525, 1526, 1527, 1528, 1529, 153, 1530, 1531, 1532, 1533, 1534, 1535, 1536,
  1537, 1538, 1539, 154, 1540, 1541, 1542, 1543, 1544, 1545, 1546, 1547, 1549,
  155, 1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558, 1559, 156, 1560,
  1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 157, 1570, 1571, 1572,
  1573, 1574, 1575, 1576, 1577, 1578, 1579, 158, 1580, 1581, 1582, 1583, 1584,
  1585, 1586, 1587, 1588, 1589, 159, 1590, 1591, 1592, 1593, 1594, 1595, 1596,
  1598, 1599, 160, 1600, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608, 1609,
  161, 1610, 1611, 1612, 1613, 1614, 1615, 1616, 1617, 1618, 1619, 162, 1620,
  1621, 1622, 1623, 1624, 1625, 1626, 1627, 1628, 1629, 163, 1630, 1631, 1632,
  1633, 1634, 1635, 1636, 1637, 1638, 1639, 164, 1640, 1641, 1642, 1643, 1644,
  1645, 1646, 1647, 1648, 1649, 165, 1650, 1651, 1652, 1653, 1654, 1655, 1656,
  1657, 1658, 1659, 166, 1660, 1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668,
  1669, 167, 1670, 1671, 1672, 1673, 1674, 1675, 1676, 1677, 1678, 1679, 168,
  1681, 1682, 1683, 1684, 1685, 1686, 1687, 1688, 1689, 169, 1690, 1691, 1692,
  1693, 1694, 1695, 1696, 1697, 1698, 1699, 170, 1700, 1701, 1702, 1703, 1704,
  1705, 1706, 1707, 1708, 1709, 171, 1710, 1711, 1712, 1713, 1714, 1715, 1716,
  1717, 1718, 1719, 172, 1720, 1721, 1722, 1723, 1724, 1725, 1726, 1727, 1728,
  1729, 173, 1730, 1731, 1732, 1733, 1734, 1735, 1736, 1737, 1738, 1739, 174,
  1740, 1741, 1742, 1743, 1744, 1745, 1746, 1747, 1748, 1749, 175, 1750, 1751,
  1752, 1753, 1754, 1755, 1756, 1757, 1758, 1759, 176, 1760, 1761, 1762, 1763,
  1764, 1765, 1766, 1767, 1768, 1769, 177, 1770, 1771, 1772, 1773, 1774, 1775,
  1776, 1777, 1778, 1779, 178, 1780, 1781, 1782, 1783, 1784, 1785, 1786, 1787,
  1788, 1789, 179, 1790, 1791, 1792, 1793, 1794, 1795, 1796, 1797, 1798, 1799,
  180, 1800, 1801, 1802, 1803, 1804, 1805, 1806, 1807, 1808, 1809, 181, 1810,
  1811, 1812, 1813, 1814, 1815, 1816, 1817, 1818, 1819, 182, 1820, 1821, 1822,
  1823, 1824, 1825, 1826, 1827, 1828, 1829, 183, 1830, 1831, 1832, 1833, 1834,
  1835, 1836, 1837, 1838, 1839, 184, 1840, 1841, 1842, 1843, 1844, 1845, 1846,
  1847, 1848, 1849, 185, 1850, 1851, 1852, 1853, 1854, 1855, 1856, 1857, 1858,
  1859, 186, 1860, 1861, 1863, 1864, 1865, 1866, 1867, 1868, 1869, 187, 1870,
  1871, 1872, 1873, 1874, 1875, 1876, 1877, 1878, 1879, 188, 1880, 1881, 1882,
  1883, 1884, 1885, 1886, 1887, 1888, 1889, 189, 1890, 1891, 1892, 1893, 1894,
  1895, 1896, 1897, 1898, 1899, 190, 1900, 1901, 1902, 1903, 1904, 1905, 1906,
  1907, 1908, 1909, 191, 1910, 1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918,
  1919, 192, 1920, 1921, 1922, 1923, 1924, 1925, 1926, 1927, 1928, 1929, 193,
  1930, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938, 1939, 194, 1940, 1941,
  1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 195, 1950, 1951, 1952, 1953,
  1955, 1956, 1957, 1958, 1959, 196, 1960, 1961, 1962, 1963, 1964, 1965, 1966,
  1967, 1968, 1969, 197, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978,
  1979, 198, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 199,
  1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 200, 2000, 2001,
  2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 201, 2010, 2011, 2012, 2013,
  2014, 2015, 2016, 2017, 2018, 2019, 202, 2020, 2021, 2022, 2023, 2024, 2025,
  2026, 2027, 2028, 2029, 203, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037,
  2038, 204, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 205,
  2050, 2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 206, 2060, 2061,
  2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 207, 2070, 2071, 2072, 2073,
  2074, 2075, 2076, 2077, 2078, 2079, 208, 2080, 2081, 2082, 2083, 2084, 2085,
  2086, 2087, 2088, 2089, 209, 2090, 2091, 2092, 2093, 2094, 2095, 2096, 2097,
  2098, 2099, 210, 2100, 2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109,
  211, 2110, 2111, 2112, 2113, 2114, 2115, 2116, 2117, 2118, 2119, 212, 2120,
  2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128, 2129, 213, 2130, 2131, 2132,
  2133, 2134, 2135, 2136, 2137, 2138, 2139, 214, 2140, 2141, 2142, 2143, 2144,
  2145, 2146, 2147, 2148, 2149, 215, 2150, 2151, 2152, 2153, 2154, 2155, 2156,
  2157, 2158, 2159, 216, 2160, 2161, 2162, 2163, 2164, 2165, 2166, 2167, 2168,
  2169, 2170, 2172, 2174, 2175, 2176, 2177, 2178, 2179, 218, 2180, 2181, 2182,
  2183, 2184, 2185, 2186, 2187, 2188, 2189, 219, 2190, 2191, 2192, 2193, 2194,
  2195, 2196, 2197, 2198, 2199, 220, 2200, 2201, 2202, 2203, 2204, 2205, 2206,
  2207, 2208, 2209, 221, 2210, 2212, 2213, 2214, 2215, 2216, 2217, 2218, 222,
  2220, 2221, 2222, 2223, 2224, 2225, 2226, 2227, 2228, 2229, 223, 2230, 2231,
  2232, 2233, 2234, 2235, 2236, 2237, 2238, 2239, 224, 2240, 2241, 2243, 2244,
  2245, 2246, 2247, 2249, 225, 2250, 2251, 2252, 2253, 2254, 2255, 2256, 2257,
  2258, 226, 2260, 2261, 2263, 2264, 2265, 2267, 2268, 2269, 227, 2270, 2271,
  2272, 2273, 2274, 2275, 2276, 2277, 2278, 2279, 228, 2280, 2281, 2282, 2283,
  2284, 2285, 2286, 2287, 2288, 2289, 229, 2290, 2291, 2292, 2293, 2294, 2295,
  2296, 2298, 2299, 230, 2301, 2302, 2303, 2304, 2306, 2307, 2308, 2309, 231,
  2310, 2311, 2312, 2313, 2314, 2315, 2316, 2317, 2318, 2319, 232, 2320, 2321,
  2322, 2323, 2324, 2325, 2326, 2327, 2328, 2329, 233, 2330, 2331, 2332, 2333,
  2334, 2335, 2336, 2337, 2338, 2339, 234, 2340, 2341, 2342, 2343, 2344, 2345,
  2346, 2347, 2348, 2349, 235, 2350, 2351, 2352, 2354, 2355, 2356, 2357, 2358,
  2359, 236, 2360, 2361, 2362, 2363, 2364, 2365, 2367, 2368, 237, 2370, 2371,
  2372, 2373, 2374, 2375, 2378, 2379, 238, 2380, 2381, 2382, 2383, 2384, 2385,
  2386, 2387, 2388, 2389, 239, 2390, 2391, 2392, 2393, 2394, 2395, 2396, 2397,
  2398, 2399, 240, 2400, 2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409,
  241, 2410, 2411, 2413, 2414, 2415, 2416, 2417, 2418, 2419, 242, 2420, 2421,
  2422, 2423, 2424, 2425, 2426, 2427, 2428, 2429, 243, 2430, 2431, 2432, 2433,
  2434, 2435, 2436, 2437, 2438, 2439, 244, 2440, 2441, 2442, 2443, 2444, 2445,
  2446, 2447, 2448, 2449, 245, 2450, 2451, 2452, 2453, 2454, 2455, 2456, 2458,
  2459, 246, 2460, 2461, 2462, 2463, 2464, 2465, 2466, 2467, 2468, 2469, 247,
  2470, 2471, 2472, 2473, 2474, 2475, 2476, 2477, 2478, 2479, 248, 2480, 2481,
  2482, 2483, 2484, 2485, 2486, 2487, 2488, 2489, 249, 2490, 2491, 2492, 2493,
  2494, 2495, 2496, 2497, 2498, 2499, 250, 2500, 251, 252, 253, 254, 255, 256,
  257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
  273, 274, 275, 276, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289,
  290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305,
  306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321,
  322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337,
  338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353,
  354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369,
  370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385,
  386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401,
  402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417,
  418, 419, 420, 421, 422, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434,
  435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450,
  451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 463, 464, 465, 466, 467,
  468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483,
  484, 485, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500,
  501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
  517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532,
  533, 534, 535, 536, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549,
  550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565,
  566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581,
  582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597,
  599, 600, 601, 602, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615,
  616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631,
  632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647,
  648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663,
  664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679,
  680, 681, 682, 683, 684, 685, 686, 688, 689, 690, 691, 693, 694, 695, 696, 697,
  698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713,
  714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729,
  730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745,
  746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761,
  762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777,
  778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793,
  794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809,
  810, 811, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 826,
  827, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842,
  843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858,
  859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874,
  875, 876, 877, 878, 879, 880, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891,
  892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 907, 908,
  909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924,
  925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 940,
  941, 942, 943, 944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955, 956,
  957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972,
  973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988,
  989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999
};*/ //nruns = 2355
//Task configuration
TString cutFile="./cutfile/D0DsDplusCuts.root";          // file containing the cuts for the different mesons
  														 // to generate the cut file: 1) move to cutfile directory
  														 //                           2) .L makeCutsTreeCreator.C
  														 //                           3) makeCutsTreeCreator();
  														 // to run with ROOT5/6 generate the cut file using AliPhysics built on ROOT5/6


bool isOutTree=true;
bool isOutNorm=false;

//************************************

void runAnalysis_v2()
{

     // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    // Load libraries

 // MC generator libraries
    //gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGJE -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/PYTHIA6 -g");

    //gSystem->Load("liblhapdf.so");
    //gSystem->Load("libpythia6.4.25.so");
    //gSystem->Load("libEGPythia6");
    //gSystem->Load("libAliPythia6");


    //gSystem->Load("libpythia6");

    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = runLocal;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = runGridTest;

    // since we will compile a class, tell root where to look for headers
// #if !defined (__CINT__) || defined (__CLING__)
//     gInterpreter->ProcessLine(".include $ROOTSYS/include");
//     gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
// #else
//     gROOT->ProcessLine(".include $ROOTSYS/include");
//     gROOT->ProcessLine(".include $ALICE_ROOT/include");
// #endif
//
// 	gSystem->Load("libEG");
//     gSystem->Load("libEGPythia6");
//     gSystem->Load("libgeant321.so");
//     gSystem->Load("liblhapdf.so");      // Parton density functions
//     gSystem->Load("libpythia6.so");     // Pythia
//     gSystem->Load("libAliPythia6.so");

    // create the analysis manager

   //  header location
//         gInterpreter->ProcessLine(".include $ROOTSYS/include");
//         gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    if(isOutNorm) mgr->SetCommonFileName("NormalizationCounterOutputBkg.root");
    //mgr->SetCommonFileName("AnalysisResults_recoOmegac_001.root");

    AliESDEvent *esdE = new AliESDEvent();
    esdE->CreateStdContent();
    AliESDVertex *vtx = new AliESDVertex(0.,0.,100);
    vtx->SetName("VertexTracks");
    vtx->SetTitle("VertexTracks");
    esdE->SetPrimaryVertexTracks(vtx);
    AliDummyHandler *dumH = new AliDummyHandler;//  static_cast<AliDummyHandler*>mgr->GetInputEventHandler();
    dumH->SetEvent(esdE);
    mgr->SetInputEventHandler(dumH);

    //mgr->SetInputEventHandler(new AliDummyHandler());
    AliMCEventHandler *mc = new AliMCEventHandler();
    mc->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mc);

//TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));

//AliAnalysisTaskPIDResponse* PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse*>(PIDadd.Exec());


    gInterpreter->LoadMacro("./R5Detector.cxx++g");
      gInterpreter->LoadMacro("AliAnalysisTaskSEXiccTopKpipi.cxx++g");
      AliAnalysisTaskSEXiccTopKpipi *task = reinterpret_cast<AliAnalysisTaskSEXiccTopKpipi*>(gInterpreter->ExecuteMacro("AddTaskXiccTopKpipi.C"));

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(0);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {

        // if you want to run locally, we need to define some input
        TChain* chainAOD = new TChain("TE");
        //TChain *chainAODfriend = new TChain("TreeK");

        // add a few files to the chain (change this so that your local files are added)
        chainAOD->Add(Form("%s/galice.root",pathToLocalAODfiles.Data()));

        //for(Int_t i=1; i<=50; i++){
        //	chainAOD->Add(Form("%s/%d/galice.root",pathToLocalAODfiles.Data(),i));
        //}

        chainAOD->Print();

        // start the analysis locally, reading the events from the tchain
        //mgr->StartAnalysis("local", chainAOD);
        mgr->SetCollectSysInfoEach(5);
        //mgr->StartAnalysis("local", chainAOD,25,226);//, 100);//,2,103);//, 10);//,1,0);//,2,0);
        TStopwatch timer;
        timer.Start();
        mgr->StartAnalysis("local", chainAOD);//, 100);//,2,103);//, 10);//,1,0);//,2,0);
        timer.Stop();
        printf("RealTime =%7.3f s, CpuTime =%7.3f s\n",timer.RealTime(),timer.CpuTime());


    } else {

        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/PYTHIA6 -I$ALICE_PHYSICS/PWGLF");
        //gSystem->Load("libEG");
        //gSystem->Load("libEGPythia6");
        //gSystem->Load("liblhapdf.so");      // Parton density functions
        //gSystem->Load("libpythia6.so");     // Pythia
        //gSystem->Load("libAliPythia6.so");
        // make sure your source files get copied to grid
        //alienHandler->SetAdditionalLibs("AliHFCutOptTreeHandler.cxx AliHFCutOptTreeHandler.h");
        if(isOutTree){
          alienHandler->SetAdditionalLibs("R5Detector.cxx R5Detector.h AliAnalysisTaskWeakDecayVertexer_mod.cxx AliAnalysisTaskWeakDecayVertexer_mod.h Al AliRecoDecayOmegac.h AliRecoDecayOmegacc.cxx AliRecoDecayOmegacc.h AliRecoDecayOmegaccc.cxx AliRecoDecayOmegaccc.h AliAnalysisTaskSEOmegacccToOmega3Pi.cxx AliAnalysisTaskSEXiccTopKpipi.cxx");
          alienHandler->SetAnalysisSource("R5Detector.cxx AliAnalysisTaskWeakDecayVertexer_mod.cxx AliRecoDecayOmegac.cxx AliRecoDecayOmegacc.cxx AliRecoDecayOmegaccc.cxx AliAnalysisTaskSEOmegacccToOmega3Pi.cxx");
        }
        if(isOutNorm){
          alienHandler->SetAdditionalLibs("R5Detector.cxx R5Detector.h AliAnalysisTaskWeakDecayVertexer_mod.cxx AliAnalysisTaskWeakDecayVertexer_mod.h AliRecoDecayOmegac.cxx AliRecoDecayOmegac.h AliRecoDecayOmegacc.cxx AliRecoDecayOmegacc.h AliRecoDecayOmegaccc.cxx AliRecoDecayOmegaccc.h AliAnalysisTaskSENormCounterOmegaccc.cxx AliAnalysisTaskSENormCounterOmegaccc.h");
          alienHandler->SetAnalysisSource("R5Detector.cxx AliAnalysisTaskWeakDecayVertexer_mod.cxx AliRecoDecayOmegac.cxx AliRecoDecayOmegacc.cxx AliRecoDecayOmegaccc.cxx AliAnalysisTaskSENormCounterOmegaccc.cxx");
        }

        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion(aliPhysVersion.Data());

        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // select the input data
        alienHandler->SetGridDataDir(gridDataDir.Data());
        //alienHandler->SetDataPattern(Form("%s/*AliAOD.root",gridDataPattern.Data()));
        //alienHandler->SetFriendChainName("AliAOD.VertexingHF.root");
        alienHandler->SetDataPattern(Form("%s*galice.root",gridDataPattern.Data()));
        alienHandler->SetTreeName("TE");
        //alienHandler->SetFriendChainName("Kinematics.root");

        // MC has no prefix, data has prefix 000
        //if(!isRunOnMC)alienHandler->SetRunPrefix("000");

        //alienHandler->AddDataFile("/alice/cern.ch/user/a/afestant/testHIJINGPYTHIAgen_omegaccc_bkg/background_analysis_Kr/input_data_Kr.xml");
        //alienHandler->SetUseMCchain();
        //alienHandler->SetNrunsPerMaster(1);
        //if(!isRunOnMC)
        //alienHandler->SetRunPrefix("0");
        for(Int_t k=0; k<nruns; k++){
            alienHandler->AddRunNumber(runlist[k]);

        }
        alienHandler->SetNrunsPerMaster(10);

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(1);
        alienHandler->SetExecutable("myTask_config3.sh");

        // specify how many seconds your job may take
        alienHandler->SetTTL(70000);
        alienHandler->SetJDLName("myTask_config3.jdl");

        //alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        if(isOutNorm) alienHandler->SetOutputFiles("AnalysisNormalizationResults.root");
        alienHandler->SetMaxMergeStages(3); //2, 3
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        alienHandler->SetGridWorkingDir(gridWorkingDir.Data());
        alienHandler->SetGridOutputDir(gridOutputDir.Data());

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if(gridTest) {

            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        }
        else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(runMode.Data()); //terminate
            mgr->StartAnalysis("grid");
        }
    }
}
