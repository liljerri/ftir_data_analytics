"""
Peak definitions based on literature values

This file provides peak definitions based on different literature sources. The
source of the definition is provides in the comments of this file, and the
name of the definition follows the same convention:
    firstauthor_solvent_year

The peak definitions are dictionary objects with four items

means: List
    list of the mean frequencies for a given peak

uncertainties: List of tuples
    list of tuples providing to lower and upper bounds of the absorbance band

relative_uncertainties : List
    list of the plus/minus values around the mean. The `relative_uncertainties`
    are used to calculate the `uncertainties`.

assignments : List
    list of the literature peak assignments

"""


# Dong et. al., Biochemistry. 1990
dong_h2o_1990 = {
    'means': [1624, 1627, 1632, 1638, 1642, 1650, 1656, 1666, 1672, 1680,
              1688],
    'uncertainties': [(1623.5, 1624.5), (1626, 1628), (1631, 1633),
                      (1637, 1639), (1641, 1643), (1649, 1651), (1654, 1658),
                      (1665, 1667), (1671, 1673), (1679, 1681), (1687, 1689)],
    'relative_uncertainties': [0.5, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1],
    'assignments': ['beta_sheet', 'beta_sheet',
                    'beta_sheet_and_extended_chain', 'beta_sheet', 'unordered',
                    'alpha_helix', 'turn', 'turn', 'turn', 'turn']
}


# Yang et. al., Nature Protocols. 2015
yang_h20_2015 = {
    'means': [1624, 1627, 1633, 1638, 1642, 1648, 1656, 1663, 1667, 1675, 1680,
              1685, 1691, 1696],
    'uncertainties': [(1623, 1625), (1625, 1629), (1631, 1635), (1636, 1640),
                      (1641, 1643), (1646, 1650), (1654, 1658), (1660, 1666),
                      (1666, 1668), (1674, 1676), (1678, 1682), (1683, 1687),
                      (1689, 1693), (1694, 1698)],
    'relative_uncertainties': [1, 2, 2, 2, 1, 2, 2, 3, 1, 1, 2, 2, 2, 2],
    'assignments': ['\u03B2-sheet', '\u03B2-sheet', '\u03B2-sheet', '\u03B2-sheet',
                    '\u03B2-sheet', 'unordered', '\u03B1-helix', '310-helix',
                    '\u03B2-turn', '\u03B2-turn', '\u03B2-turn', '\u03B2-turn',
                    '\u03B2-sheet', '\u03B2-sheet']
}


yang_d20_2015 = {
    'means': [1624, 1631, 1637, 1641, 1645, 1653, 1663, 1671, 1675, 1683, 1689,
              1694],
    'uncertainties': [(1620, 1628), (1628, 1634), (1634, 1640), (1639, 1643),
                      (1641, 1649), (1649, 1657), (1659, 1667), (1668, 1674),
                      (1670, 1680), (1681, 1685), (1687, 1691), (1692, 1696)],
    'relative_uncertainties': [4, 3, 3, 2, 4, 4, 4, 3, 5, 2, 2, 2],
    'assignments': ['beta_sheet', 'beta_sheet', 'beta_sheet', '310_helix',
                    'unordered', 'alpha_helix', 'beta_turn', 'beta_turn',
                    'beta_sheet', 'beta_turn', 'beta_turn', 'beta_turn']
}
