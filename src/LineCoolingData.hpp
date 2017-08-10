/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file LineCoolingData.hpp
 *
 * @brief Internal representation of the external data used for line cooling:
 * header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LINECOOLINGDATA_HPP
#define LINECOOLINGDATA_HPP

#include "LineCoolingDataLocation.hpp"
#include <string>

/**
 * @brief Names of supported five level elements.
 */
enum LineCoolingDataFiveLevelElement {
  /*! @brief Nitrogen I. */
  NI = 0, // done
  /*! @brief Nitrogen II. */
  NII, // done
  /*! @brief Oxygen I. */
  OI, // done
  /*! @brief Oxygen II. */
  OII, // done
  /*! @brief Oxygen III. */
  OIII, // done
  /*! @brief Neon III. */
  NeIII, // done
  /*! @brief Sulphur II. */
  SII, // done
  /*! @brief Sulphur III. */
  SIII,
  /*! @brief Carbon II. */
  CII, // done
  /*! @brief Carbon III. */
  CIII, // done
  /*! @brief Counter. Should always be the last element! */
  LINECOOLINGDATA_NUMFIVELEVELELEMENTS
};

/**
 * @brief Names of supported two level elements.
 */
enum LineCoolingDataTwoLevelElement {
  /*! @brief Nitrogen III. */
  NIII = 0,
  /*! @brief Neon II. */
  NeII,
  /*! @brief Counter. Should always be the last element! */
  LINECOOLINGDATA_NUMTWOLEVELELEMENTS
};

/**
 * @brief Names of the two level element data fields.
 */
enum LineCoolingDataTwoLevelFields {
  /*! @brief Energy difference between the ground level and first excited level
   *  (in K). */
  TWOLEVELFIELD_ENERGY_DIFFERENCE,
  /*! @brief Transition probability for deexcitation from the first excited
   *  level to the ground level (in s^-1). */
  TWOLEVELFIELD_TRANSITION_PROBABILITY,
  /*! @brief Velocity-averaged collision strength for the transition from the
   *  ground level to the first excited level (at 10,000 K). */
  TWOLEVELFIELD_COLLISION_STRENGTH,
  /*! @brief Statistical weight of the ground level. */
  TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_0,
  /*! @brief Statistical weight of the first excited level. */
  TWOLEVELFIELD_INVERSE_STATISTICAL_WEIGHT_1,
  /*! @brief Counter. Should always be the last element! */
  LINECOOLINGDATA_NUMTWOLEVELFIELDS
};

/**
 * @brief Convenient names for transitions between levels.
 */
enum LineCoolingDataTransition {
  /*! @brief A transition from the ground level to the first excited level. */
  TRANSITION_0_to_1 = 0,
  /*! @brief A transition from the ground level to the second excited level. */
  TRANSITION_0_to_2,
  /*! @brief A transition from the ground level to the third excited level. */
  TRANSITION_0_to_3,
  /*! @brief A transition from the ground level to the fourth excited level. */
  TRANSITION_0_to_4,
  /*! @brief A transition from the first excited level to the second excited
   *  level. */
  TRANSITION_1_to_2,
  /*! @brief A transition from the first excited level to the third excited
   *  level. */
  TRANSITION_1_to_3,
  /*! @brief A transition from the first excited level to the fourth excited
   *  level. */
  TRANSITION_1_to_4,
  /*! @brief A transition from the second excited level to the third excited
   *  level. */
  TRANSITION_2_to_3,
  /*! @brief A transition from the second excited level to the fourth excited
   *  level. */
  TRANSITION_2_to_4,
  /*! @brief A transition from the third excited level to the fourth excited
   *  level. */
  TRANSITION_3_to_4,
  /*! @brief Counter. Should always be the last element! */
  NUMBER_OF_TRANSITIONS
};

/**
 * @brief Internal representation of the line cooling data in "atom4.dat".
 *
 * The cooling by collisionally excited line radiation is based on section 3.5
 * of Osterbrock, D. E. & Ferland, G. J. 2006, Astrophysics of Gaseous Nebulae
 * and Active Galactic Nuclei, 2nd edition
 * (http://adsabs.harvard.edu/abs/2006agna.book.....O).
 *
 * The data used comes from Pradhan A. K. & Peng, J. 1995, The Analsis of
 * Emission Lines, STScI Symp. 8, Cambridge Univ. Press, and data that used to
 * be available from http://www-astronomy.mps.ohio-state.edu/∼pradhan/.
 * Unfortunately, these data are no longer accessible, so we just hope our
 * legacy version is accurate.
 *
 * We use data from
 *  - http://cdsweb.u-strasbg.fr/tipbase/home.html
 *  - Blum, R. D. & Pradhan, A. K. 1992, ApJS, 80, 425
 *    (http://adsabs.harvard.edu/abs/1992ApJS...80..425B)
 *  - Galavis, M. E., Mendoza, C. & Zeippen, C. J. 1998, A&AS, 131, 499
 *    (http://adsabs.harvard.edu/abs/1998A%26AS..131..499G)
 *  - Saraph, H. E. & Tully, J. A. 1994, A&AS, 107, 29
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..107...29S)
 *  - Kaufman, V. & Sugar, J. 1986, JPCRD, 15, 321
 *    (http://adsabs.harvard.edu/abs/1986JPCRD..15..321K)
 *  - Griffin, D. C, Mitnik, D. M., Badnell, N. R. 2001, JPhB, 34, 4401
 *    (http://adsabs.harvard.edu/abs/2001JPhB...34.4401G)
 *  - Galavis, M. E., Mendoza, C. & Zeippen, C. J. 1997, A&AS, 123, 159
 *    (http://adsabs.harvard.edu/abs/1997A%26AS..123..159G)
 *  - Lennon, D. J. & Burke, V. M. 1994, A&AS, 103, 273
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..103..273L)
 *  - Butler, K. & Zeippen, C. J. 1994, A&AS, 108, 1
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..108....1B)
 *  - Zatsarinny, O. & Tayal, S. S. 2003, ApJS, 148, 575
 *    (http://adsabs.harvard.edu/abs/2003ApJS..148..575Z)
 *  - Berrington, K. A. 1988, JPhB, 21, 1083
 *    (http://adsabs.harvard.edu/abs/1988JPhB...21.1083B)
 *  - Froese Fischer, C. & Tachiev, G. 2004, ADNDT, 87, 1
 *    (http://adsabs.harvard.edu/abs/2004ADNDT..87....1F)
 *  - Tayal, S. S. 2000, ADNDT, 76, 191
 *    (http://adsabs.harvard.edu/abs/2000ADNDT..76..191T)
 *  - Kisielius, R., Storey, P. J., Ferland, G. J. & Keenan, F. P. 2009, MNRAS,
 *    397, 903 (http://adsabs.harvard.edu/abs/2009MNRAS.397..903K)
 *  - Tayal, S. S. 2008, A&A, 486, 629
 *    (http://adsabs.harvard.edu/abs/2008A%26A...486..629T)
 *  - Berrington, K. A., Burke, P. G., Dufton, P. L. & Kingston, A. E. 1985,
 *    ADNDT, 33, 195 (http://adsabs.harvard.edu/abs/1985ADNDT..33..195B)
 *  - Tayal, S. S. & Zatsarinny, O. 2010, ApJS, 188, 32
 *    (http://adsabs.harvard.edu/abs/2010ApJS..188...32T)
 *  - Mendoza, C. & Zeippen, C. J. 1982, MNRAS, 199, 1025
 *    (http://adsabs.harvard.edu/abs/1982MNRAS.199.1025M)
 *  - Hudson, C. E., Ramsbottom, C. A. & Scott, M. P. 2012, ApJ, 750, 65
 *    (http://adsabs.harvard.edu/abs/2012ApJ...750...65H)
 */
class LineCoolingData {
private:
  /*! @brief Velocity-averaged collision strengths at 10,000 K. */
  double _collision_strength[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                            [NUMBER_OF_TRANSITIONS];

  /*! @brief Exponents for the temperature variation of the collision
   *  strengths. */
  double _collision_strength_exponent[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                     [NUMBER_OF_TRANSITIONS];

  /*! @brief Transition probabilities for deexcitation between different
   *  levels. */
  double _transition_probability[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                [NUMBER_OF_TRANSITIONS];

  /*! @brief Energy differences for the transitions between different levels
   *  (in K). */
  double _energy_difference[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                           [NUMBER_OF_TRANSITIONS];

  /*! @brief Inverse statistical weights for the different levels. */
  double _inverse_statistical_weight[LINECOOLINGDATA_NUMFIVELEVELELEMENTS][5];

  /*! @brief Data values for the two level elements. */
  double _two_level_element_data[LINECOOLINGDATA_NUMTWOLEVELELEMENTS]
                                [LINECOOLINGDATA_NUMTWOLEVELFIELDS];

  static bool read_values(std::string line, double *array, unsigned int size);

public:
  LineCoolingData();

  double get_collision_strength(LineCoolingDataFiveLevelElement element,
                                LineCoolingDataTransition transition) const;
  double
  get_collision_strength_exponent(LineCoolingDataFiveLevelElement element,
                                  LineCoolingDataTransition transition) const;
  double get_transition_probability(LineCoolingDataFiveLevelElement element,
                                    LineCoolingDataTransition transition) const;
  double get_energy_difference(LineCoolingDataFiveLevelElement element,
                               LineCoolingDataTransition transition) const;
  double get_statistical_weight(LineCoolingDataFiveLevelElement element,
                                unsigned char level) const;

  static int simq(double A[5][5], double B[5]);

  double get_cooling(double temperature, double electron_density,
                     const double *abundances) const;

  void linestr(double temperature, double electron_density,
               const double *abundances, double &c6300, double &c9405,
               double &c6312, double &c33mu, double &c19mu, double &c3729,
               double &c3727, double &c7330, double &c4363, double &c5007,
               double &c52mu, double &c88mu, double &c5755, double &c6584,
               double &c4072, double &c6717, double &c6725, double &c3869,
               double &cniii57, double &cneii12, double &cneiii15,
               double &cnii122, double &cii2325, double &ciii1908,
               double &coii7325, double &csiv10) const;
};

#endif // LINECOOLINGDATA_HPP
