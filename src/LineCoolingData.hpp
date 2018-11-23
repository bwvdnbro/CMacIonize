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

#include <string>
#include <vector>

/**
 * @brief Names of supported five level elements.
 *
 * The order is historical: this was the order of the elements in Kenny's code.
 */
enum LineCoolingDataFiveLevelElement {
  /*! @brief Nitrogen I. */
  NI = 0,
  /*! @brief Nitrogen II. */
  NII,
  /*! @brief Oxygen I. */
  OI,
  /*! @brief Oxygen II. */
  OII,
  /*! @brief Oxygen III. */
  OIII,
  /*! @brief Neon III. */
  NeIII,
  /*! @brief Sulphur II. */
  SII,
  /*! @brief Sulphur III. */
  SIII,
  /*! @brief Carbon II. */
  CII,
  /*! @brief Carbon III. */
  CIII,
  /*! @brief Counter. Should always be the last element! */
  LINECOOLINGDATA_NUMFIVELEVELELEMENTS
};

/**
 * @brief Names of supported two level elements.
 *
 * Note that we start counting from the number of five level elements!
 */
enum LineCoolingDataTwoLevelElement {
  /*! @brief Nitrogen III. */
  NIII = LINECOOLINGDATA_NUMFIVELEVELELEMENTS,
  /*! @brief Neon II. */
  NeII,
  /*! @brief Sulphur IV. */
  SIV,
  /*! @brief Counter. Should always be the last element! */
  LINECOOLINGDATA_NUMELEMENTS
};

/*! @brief Number of two level elements, which is defined as the total number of
 *  elements minus the number of five level elements. */
#define LINECOOLINGDATA_NUMTWOLEVELELEMENTS                                    \
  (LINECOOLINGDATA_NUMELEMENTS - LINECOOLINGDATA_NUMFIVELEVELELEMENTS)

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
 * We use data from a large number of different sources, many of which are part
 * of the the IRON project and were located using their online database
 * (http://cdsweb.u-strasbg.fr/tipbase/home.html). We also used the extensive
 * list of data sources in Lykins, M. L., Ferland, G. J., Kisielius, R.,
 * Chatzikos, M., Porter, R. L., van Hoof, P. A. M., Williams, R. J. R., Keenan,
 * F. P. & Stancil, P. C. 2015, ApJ, 807, 118
 * (http://adsabs.harvard.edu/abs/2015ApJ...807..118L) to locate data for ions
 * that seem to be missing from the IRON database.
 *
 * The actual data used comes from (in alphabetical order):
 *  - Berrington, K. A. 1988, JPhB, 21, 1083
 *    (http://adsabs.harvard.edu/abs/1988JPhB...21.1083B) (OI)
 *  - Berrington, K. A., Burke, P. G., Dufton, P. L. & Kingston, A. E. 1985,
 *    ADNDT, 33, 195 (http://adsabs.harvard.edu/abs/1985ADNDT..33..195B) (CIII)
 *  - Blum, R. D. & Pradhan, A. K. 1992, ApJS, 80, 425
 *    (http://adsabs.harvard.edu/abs/1992ApJS...80..425B) (NIII)
 *  - Butler, K. & Zeippen, C. J. 1994, A&AS, 108, 1
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..108....1B) (NeIII)
 *  - Froese Fischer, C. & Tachiev, G. 2004, ADNDT, 87, 1
 *    (http://adsabs.harvard.edu/abs/2004ADNDT..87....1F) (NI, OII, CII, CIII)
 *  - Galavis, M. E., Mendoza, C. & Zeippen, C. J. 1997, A&AS, 123, 159
 *    (http://adsabs.harvard.edu/abs/1997A%26AS..123..159G) (NII, OI, OIII,
 *    NeIII)
 *  - Galavis, M. E., Mendoza, C. & Zeippen, C. J. 1998, A&AS, 131, 499
 *    (http://adsabs.harvard.edu/abs/1998A%26AS..131..499G) (NIII)
 *  - Griffin, D. C, Mitnik, D. M., Badnell, N. R. 2001, JPhB, 34, 4401
 *    (http://adsabs.harvard.edu/abs/2001JPhB...34.4401G) (NeII)
 *  - Hudson, C. E., Ramsbottom, C. A. & Scott, M. P. 2012, ApJ, 750, 65
 *    (http://adsabs.harvard.edu/abs/2012ApJ...750...65H) (SIII)
 *  - Kaufman, V. & Sugar, J. 1986, JPCRD, 15, 321
 *    (http://adsabs.harvard.edu/abs/1986JPCRD..15..321K) (NeII)
 *  - Kisielius, R., Storey, P. J., Ferland, G. J. & Keenan, F. P. 2009, MNRAS,
 *    397, 903 (http://adsabs.harvard.edu/abs/2009MNRAS.397..903K) (OII)
 *  - Lennon, D. J. & Burke, V. M. 1994, A&AS, 103, 273
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..103..273L) (NII, OIII)
 *  - Martin, W. C., Zalubas, R. & Musgrove, A. 1990, JPCRD, 19, 821
 *    (http://adsabs.harvard.edu/abs/1990JPCRD..19..821M) (SIV)
 *  - Mendoza, C. & Zeippen, C. J. 1982, MNRAS, 199, 1025
 *    (http://adsabs.harvard.edu/abs/1982MNRAS.199.1025M) (SIII)
 *  - Pradhan, A. 1995, The Analysis of Emission Lines: A Meeting in Honor of
 *    the 70th Birthdays of D. E. Osterbrock & M. J. Seaton, proceedings of the
 *    Space Telescope Science Institute Symposium, held in Baltimore, Maryland
 *    May 16--18, 1994, Eds.: Robert Williams and Mario Livio, Cambridge
 *    University Press, p. 8.
 *    (http://adsabs.harvard.edu/abs/1995aelm.conf....8P) (SIV)
 *  - Saraph, H. E. & Storey, P. J. 1999, A&AS, 134, 369
 *    (http://adsabs.harvard.edu/abs/1999A%26AS..134..369S) (SIV)
 *  - Saraph, H. E. & Tully, J. A. 1994, A&AS, 107, 29
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..107...29S) (NeII)
 *  - Tayal, S. S. 2000, ADNDT, 76, 191
 *    (http://adsabs.harvard.edu/abs/2000ADNDT..76..191T) (NI)
 *  - Tayal, S. S. 2008, A&A, 486, 629
 *    (http://adsabs.harvard.edu/abs/2008A%26A...486..629T) (CII)
 *  - Tayal, S. S. & Zatsarinny, O. 2010, ApJS, 188, 32
 *    (http://adsabs.harvard.edu/abs/2010ApJS..188...32T) (SII)
 *  - Zatsarinny, O. & Tayal, S. S. 2003, ApJS, 148, 575
 *    (http://adsabs.harvard.edu/abs/2003ApJS..148..575Z) (OI)
 *
 * Adding new lines is straightforward: add an entry in the corresponding enum
 * (LineCoolingDataFiveLevelElement or LineCoolingDataTwoLevelElement; before
 * the counter element), and initialize the data in the constructor. To compute
 * line strengths, add relevant code to linestr(). Adding new elements will
 * break some unit tests in testLineCoolingData, but should work fine.
 */
class LineCoolingData {
private:
  /*! @brief Collision strength fit parameters for the five level elements. */
  double _five_level_collision_strength[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                       [NUMBER_OF_TRANSITIONS][7];

  /*! @brief Transition probabilities for deexcitation between different
   *  levels for the five level elements. */
  double
      _five_level_transition_probability[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                        [NUMBER_OF_TRANSITIONS];

  /*! @brief Energy differences for the transitions between different levels
   *  for the five level elements (in K). */
  double _five_level_energy_difference[LINECOOLINGDATA_NUMFIVELEVELELEMENTS]
                                      [NUMBER_OF_TRANSITIONS];

  /*! @brief Inverse statistical weights for the different levels for the five
   *  level elements. */
  double _five_level_inverse_statistical_weight
      [LINECOOLINGDATA_NUMFIVELEVELELEMENTS][5];

  /*! @brief Collision strength fit parameters for the two level elements. */
  double _two_level_collision_strength[LINECOOLINGDATA_NUMTWOLEVELELEMENTS][7];

  /*! @brief Transition probabilities for deexcitation between different
   *  levels for the two level elements. */
  double _two_level_transition_probability[LINECOOLINGDATA_NUMTWOLEVELELEMENTS];

  /*! @brief Energy differences for the transitions between different levels
   *  for the two level elements (in K). */
  double _two_level_energy_difference[LINECOOLINGDATA_NUMTWOLEVELELEMENTS];

  /*! @brief Inverse statistical weights for the different levels for the two
   *  level elements. */
  double
      _two_level_inverse_statistical_weight[LINECOOLINGDATA_NUMTWOLEVELELEMENTS]
                                           [2];

  /*! @brief Prefactor for collision strengths: \f$\frac{h^2}{\sqrt{k}
   *  \left(2\pi{}m_e\right)^\frac{3}{2}\f$ (in K^0.5 m^3 s^-1). */
  double _collision_strength_prefactor;

  void compute_level_populations(int_fast32_t element,
                                 double collision_strength_prefactor, double T,
                                 double Tinv, double logT,
                                 double level_populations[5]) const;

  double compute_level_population(int_fast32_t element,
                                  double collision_strength_prefactor, double T,
                                  double Tinv, double logT) const;

public:
  LineCoolingData();

  double get_transition_probability(int_fast32_t element,
                                    int_fast32_t transition) const;
  double get_energy_difference(int_fast32_t element,
                               int_fast32_t transition) const;
  double get_statistical_weight(int_fast32_t element, uint_fast8_t level) const;

  static int solve_system_of_linear_equations(double A[5][5], double B[5]);

  double
  get_cooling(double temperature, double electron_density,
              const double abundances[LINECOOLINGDATA_NUMELEMENTS]) const;

  std::vector< std::vector< double > > get_line_strengths(
      double temperature, double electron_density,
      const double abundances[LINECOOLINGDATA_NUMELEMENTS]) const;
};

#endif // LINECOOLINGDATA_HPP
