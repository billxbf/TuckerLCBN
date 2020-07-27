/*
 * driver.cpp
 *
 *  Created on: Jun 3, 2016
 *      Author: Alicia Klinvex (amklinv@sandia.gov)
 */

#include "Tucker.hpp"
#include <iostream>
#include <cstdlib>

int main()
{
  Tucker::SizeArray* size =
      Tucker::MemoryManager::safe_new<Tucker::SizeArray>(3);
  (*size)[0] = 2;
  (*size)[1] = 3;
  (*size)[2] = 5;
  Tucker::Tensor* t =
      Tucker::MemoryManager::safe_new<Tucker::Tensor>(*size);
  Tucker::MemoryManager::safe_delete<Tucker::SizeArray>(size);
  double* data = t->data();
  data[0] = 2;    data[1] = 3;    data[2] = 5;    data[3] = 7;
  data[4] = 11;   data[5] = 13;   data[6] = 17;   data[7] = 19;
  data[8] = 23;   data[9] = 29;   data[10] = 31;  data[11] = 37;
  data[12] = 41;  data[13] = 43;  data[14] = 47;  data[15] = 53;
  data[16] = 59;  data[17] = 61;  data[18] = 67;  data[19] = 71;
  data[20] = 73;  data[21] = 79;  data[22] = 83;  data[23] = 97;
  data[24] = 101; data[25] = 103; data[26] = 107; data[27] = 109;
  data[28] = 113; data[29] = 127;

  Tucker::Matrix* mat =
      Tucker::MemoryManager::safe_new<Tucker::Matrix>(7,2);
  data = mat->data();
  data[0] = 131;  data[1] = 137; data[2] = 139;  data[3] = 149;
  data[4] = 151;  data[5] = 157; data[6] = 163;  data[7] = 167;
  data[8] = 173;  data[9] = 179; data[10] = 181; data[11] = 191;
  data[12] = 193; data[13] = 197;

  Tucker::Tensor* result = Tucker::ttm(t,0,mat,false);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat);
  data = result->data();
  double* trueData = Tucker::MemoryManager::safe_new_array<double>(105);
  trueData[0] = 763; trueData[1] = 793; trueData[2] = 815; trueData[3] = 841;
  trueData[4] = 875; trueData[5] = 893; trueData[6] = 917; trueData[7] = 1824;
  trueData[8] = 1896; trueData[9] = 1948; trueData[10] = 2012; trueData[11] = 2092;
  trueData[12] = 2136; trueData[13] = 2194; trueData[14] = 3612; trueData[15] = 3756;
  trueData[16] = 3856; trueData[17] = 3992; trueData[18] = 4144; trueData[19] = 4236;
  trueData[20] = 4354; trueData[21] = 5400; trueData[22] = 5616; trueData[23] = 5764;
  trueData[24] = 5972; trueData[25] = 6196; trueData[26] = 6336; trueData[27] = 6514;
  trueData[28] = 7856; trueData[29] = 8168; trueData[30] = 8388; trueData[31] = 8676;
  trueData[32] = 9012; trueData[33] = 9208; trueData[34] = 9462; trueData[35] = 10240;
  trueData[36] = 10648; trueData[37] = 10932; trueData[38] = 11316; trueData[39] = 11748;
  trueData[40] = 12008; trueData[41] = 12342; trueData[42] = 12552; trueData[43] = 13056;
  trueData[44] = 13396; trueData[45] = 13892; trueData[46] = 14404; trueData[47] = 14736;
  trueData[48] = 15154; trueData[49] = 15008; trueData[50] = 15608; trueData[51] = 16020;
  trueData[52] = 16596; trueData[53] = 17220; trueData[54] = 17608; trueData[55] = 18102;
  trueData[56] = 17916; trueData[57] = 18636; trueData[58] = 19120; trueData[59] = 19832;
  trueData[60] = 20560; trueData[61] = 21036; trueData[62] = 21634; trueData[63] = 20634;
  trueData[64] = 21462; trueData[65] = 22022; trueData[66] = 22834; trueData[67] = 23678;
  trueData[68] = 24222; trueData[69] = 24908; trueData[70] = 22756; trueData[71] = 23668;
  trueData[72] = 24288; trueData[73] = 25176; trueData[74] = 26112; trueData[75] = 26708;
  trueData[76] = 27462; trueData[77] = 27072; trueData[78] = 28152; trueData[79] = 28900;
  trueData[80] = 29924; trueData[81] = 31060; trueData[82] = 31752; trueData[83] = 32638;
  trueData[84] = 30432; trueData[85] = 31656; trueData[86] = 32476; trueData[87] = 33692;
  trueData[88] = 34924; trueData[89] = 35736; trueData[90] = 36754; trueData[91] = 32220;
  trueData[92] = 33516; trueData[93] = 34384; trueData[94] = 35672; trueData[95] = 36976;
  trueData[96] = 37836; trueData[97] = 38914; trueData[98] = 36012; trueData[99] = 37452;
  trueData[100] = 38440; trueData[101] = 39824; trueData[102] = 41320; trueData[103] = 42252;
  trueData[104] = 43438;
  for(int i=0; i<105; i++) {
    if(data[i] != trueData[i]) {
      std::cerr << "Mode 0 failure\n";
      return EXIT_FAILURE;
    }
  }
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(result);

  Tucker::Matrix* matt =
      Tucker::MemoryManager::safe_new<Tucker::Matrix>(2,7);
  data = matt->data();
  data[0] = 131; data[1] = 137; data[2] = 139; data[3] = 149;
  data[4] = 151; data[5] = 157; data[6] = 163; data[7] = 167;
  data[8] = 173; data[9] = 179; data[10] = 181; data[11] = 191;
  data[12] = 193; data[13] = 197;

  result = Tucker::ttm(t,0,matt,true);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(matt);
  data = result->data();
  trueData[0] = 673; trueData[1] = 725; trueData[2] = 773;
  trueData[3] = 827; trueData[4] = 883; trueData[5] = 935;
  trueData[6] = 977; trueData[7] = 1614; trueData[8] = 1738;
  trueData[9] = 1854; trueData[10] = 1984; trueData[11] = 2118;
  trueData[12] = 2242; trueData[13] = 2344; trueData[14] = 3222;
  trueData[15] = 3466; trueData[16] = 3702; trueData[17] = 3964;
  trueData[18] = 4230; trueData[19] = 4474; trueData[20] = 4684;
  trueData[21] = 4830; trueData[22] = 5194; trueData[23] = 5550;
  trueData[24] = 5944; trueData[25] = 6342; trueData[26] = 6706;
  trueData[27] = 7024; trueData[28] = 6986; trueData[29] = 7518;
  trueData[30] = 8026; trueData[31] = 8592; trueData[32] = 9170;
  trueData[33] = 9702; trueData[34] = 10152; trueData[35] = 9130;
  trueData[36] = 9822; trueData[37] = 10490; trueData[38] = 11232;
  trueData[39] = 11986; trueData[40] = 12678; trueData[41] = 13272;
  trueData[42] = 11262; trueData[43] = 12106; trueData[44] = 12942;
  trueData[45] = 13864; trueData[46] = 14790; trueData[47] = 15634;
  trueData[48] = 16384; trueData[49] = 13418; trueData[50] = 14430;
  trueData[51] = 15418; trueData[52] = 16512; trueData[53] = 17618;
  trueData[54] = 18630; trueData[55] = 19512; trueData[56] = 16086;
  trueData[57] = 17290; trueData[58] = 18486; trueData[59] = 19804;
  trueData[60] = 21126; trueData[61] = 22330; trueData[62] = 23404;
  trueData[63] = 18504; trueData[64] = 19892; trueData[65] = 21264;
  trueData[66] = 22778; trueData[67] = 24300; trueData[68] = 25688;
  trueData[69] = 26918; trueData[70] = 20386; trueData[71] = 21918;
  trueData[72] = 23426; trueData[73] = 25092; trueData[74] = 26770;
  trueData[75] = 28302; trueData[76] = 29652; trueData[77] = 24162;
  trueData[78] = 25990; trueData[79] = 27762; trueData[80] = 29728;
  trueData[81] = 31722; trueData[82] = 33550; trueData[83] = 35128;
  trueData[84] = 27342; trueData[85] = 29386; trueData[86] = 31422;
  trueData[87] = 33664; trueData[88] = 35910; trueData[89] = 37954;
  trueData[90] = 39784; trueData[91] = 28950; trueData[92] = 31114;
  trueData[93] = 33270; trueData[94] = 35644; trueData[95] = 38022;
  trueData[96] = 40186; trueData[97] = 42124; trueData[98] = 32202;
  trueData[99] = 34630; trueData[100] = 37002; trueData[101] = 39628;
  trueData[102] = 42282; trueData[103] = 44710; trueData[104] = 46828;
  for(int i=0; i<105; i++) {
    if(data[i] != trueData[i]) {
      std::cerr << "Mode 0 transpose failure\n";
      return EXIT_FAILURE;
    }
  }
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(result);

  Tucker::Matrix* mat1 =
      Tucker::MemoryManager::safe_new<Tucker::Matrix>(7,3);
  data = mat1->data();
  data[0] = 131; data[1] = 137; data[2] = 139; data[3] = 149;
  data[4] = 151; data[5] = 157; data[6] = 163; data[7] = 167;
  data[8] = 173; data[9] = 179; data[10] = 181; data[11] = 191;
  data[12] = 193; data[13] = 197; data[14] = 199; data[15] = 211;
  data[16] = 223; data[17] = 227; data[18] = 229; data[19] = 233;
  data[20] = 239;

  result = Tucker::ttm(t,1,mat1,false);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat1);
  data = result->data();
  trueData[0] = 3286; trueData[1] = 4149; trueData[2] = 3460;
  trueData[3] = 4365; trueData[4] = 3626; trueData[5] = 4569;
  trueData[6] = 3700; trueData[7] = 4665; trueData[8] = 3776;
  trueData[9] = 4767; trueData[10] = 3842; trueData[11] = 4851;
  trueData[12] = 3940; trueData[13] = 4975; trueData[14] = 12237;
  trueData[15] = 14695; trueData[16] = 12849; trueData[17] = 15427;
  trueData[18] = 13393; trueData[19] = 16083; trueData[20] = 13733;
  trueData[21] = 16479; trueData[22] = 14059; trueData[23] = 16881;
  trueData[24] = 14331; trueData[25] = 17201; trueData[26] = 14711;
  trueData[27] = 17653; trueData[28] = 24961; trueData[29] = 26623;
  trueData[30] = 26197; trueData[31] = 27931; trueData[32] = 27269;
  trueData[33] = 29067; trueData[34] = 28009; trueData[35] = 29847;
  trueData[36] = 28679; trueData[37] = 30585; trueData[38] = 29255;
  trueData[39] = 31193; trueData[40] = 30043; trueData[41] = 32029;
  trueData[42] = 37485; trueData[43] = 41797; trueData[44] = 39321;
  trueData[45] = 43861; trueData[46] = 40889; trueData[47] = 45641;
  trueData[48] = 42037; trueData[49] = 46897; trueData[50] = 43067;
  trueData[51] = 48023; trueData[52] = 43947; trueData[53] = 48995;
  trueData[54] = 45139; trueData[55] = 50319; trueData[56] = 53587;
  trueData[57] = 56969; trueData[58] = 56191; trueData[59] = 59765;
  trueData[60] = 58391; trueData[61] = 62149; trueData[62] = 60067;
  trueData[63] = 63905; trueData[64] = 61565; trueData[65] = 65455;
  trueData[66] = 62837; trueData[67] = 66799; trueData[68] = 64549;
  trueData[69] = 68615;
  for(int i=0; i<70; i++) {
    if(data[i] != trueData[i]) {
      std::cerr << "Mode 1 failure: data[" << i << "] (" << data[i]
                << ") != trueData[" << i << "] (" << trueData[i]
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(result);

  Tucker::Matrix* mat1t =
      Tucker::MemoryManager::safe_new<Tucker::Matrix>(3,7);
  data = mat1t->data();
  data[0] = 131; data[1] = 137; data[2] = 139; data[3] = 149;
  data[4] = 151; data[5] = 157; data[6] = 163; data[7] = 167;
  data[8] = 173; data[9] = 179; data[10] = 181; data[11] = 191;
  data[12] = 193; data[13] = 197; data[14] = 199; data[15] = 211;
  data[16] = 223; data[17] = 227; data[18] = 229; data[19] = 233;
  data[20] = 239;

  result = Tucker::ttm(t,1,mat1t,true);
  Tucker::MemoryManager::safe_delete<Tucker::Matrix>(mat1t);
  data = result->data();
  trueData[0] = 2476; trueData[1] = 3159; trueData[2] = 2780;
  trueData[3] = 3545; trueData[4] = 3064; trueData[5] = 3907;
  trueData[6] = 3364; trueData[7] = 4287; trueData[8] = 3560;
  trueData[9] = 4545; trueData[10] = 4034; trueData[11] = 5145;
  trueData[12] = 4252; trueData[13] = 5425; trueData[14] = 9687;
  trueData[15] = 11605; trueData[16] = 10873; trueData[17] = 13019;
  trueData[18] = 11975; trueData[19] = 14341; trueData[20] = 13127;
  trueData[21] = 15717; trueData[22] = 13981; trueData[23] = 16743;
  trueData[24] = 15753; trueData[25] = 18875; trueData[26] = 16661;
  trueData[27] = 19951; trueData[28] = 20011; trueData[29] = 21373;
  trueData[30] = 22469; trueData[31] = 23987; trueData[32] = 24739;
  trueData[33] = 26413; trueData[34] = 27115; trueData[35] = 28941;
  trueData[36] = 28913; trueData[37] = 30879; trueData[38] = 32525;
  trueData[39] = 34739; trueData[40] = 34441; trueData[41] = 36775;
  trueData[42] = 30315; trueData[43] = 33607; trueData[44] = 34037;
  trueData[45] = 37737; trueData[46] = 37471; trueData[47] = 41547;
  trueData[48] = 41059; trueData[49] = 45535; trueData[50] = 43829;
  trueData[51] = 48569; trueData[52] = 49257; trueData[53] = 54617;
  trueData[54] = 52189; trueData[55] = 57849; trueData[56] = 43597;
  trueData[57] = 46079; trueData[58] = 48947; trueData[59] = 51745;
  trueData[60] = 53881; trueData[61] = 56963; trueData[62] = 59029;
  trueData[63] = 62423; trueData[64] = 63059; trueData[65] = 66625;
  trueData[66] = 70823; trueData[67] = 74869; trueData[68] = 75067;
  trueData[69] = 79337;
  for(int i=0; i<70; i++) {
    if(data[i] != trueData[i]) {
      std::cerr << "Mode 1 transpose failure\n";
      return EXIT_FAILURE;
    }
  }
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(result);

  Tucker::MemoryManager::safe_delete_array<double>(trueData,105);
  Tucker::MemoryManager::safe_delete<Tucker::Tensor>(t);

  if(Tucker::MemoryManager::curMemUsage > 0) {
    Tucker::MemoryManager::printCurrentMemUsage();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

