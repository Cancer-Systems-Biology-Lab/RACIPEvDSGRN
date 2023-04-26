function dydt = GRHL2_47To56(t, y, p)
Inh_of_GRHL2ToZEB = p(1);
Num_of_GRHL2ToZEB = p(2);
Trd_of_GRHL2ToZEB = p(3);
Inh_of_ZEBToGRHL2 = p(4);
Num_of_ZEBToGRHL2 = p(5);
Trd_of_ZEBToGRHL2 = p(6);
Inh_of_miR200ToZEB = p(7);
Num_of_miR200ToZEB = p(8);
Trd_of_miR200ToZEB = p(9);
Inh_of_ZEBTomiR200 = p(10);
Num_of_ZEBTomiR200 = p(11);
Trd_of_ZEBTomiR200 = p(12);
Inh_of_SLUGTomiR200 = p(13);
Num_of_SLUGTomiR200 = p(14);
Trd_of_SLUGTomiR200 = p(15);
Act_of_SLUGToZEB = p(16);
Num_of_SLUGToZEB = p(17);
Trd_of_SLUGToZEB = p(18);
Act_of_ZEBToZEB = p(19);
Num_of_ZEBToZEB = p(20);
Trd_of_ZEBToZEB = p(21);
Act_of_SLUGToSLUG = p(22);
Num_of_SLUGToSLUG = p(23);
Trd_of_SLUGToSLUG = p(24);
Prod_of_GRHL2 = p(25);
Deg_of_GRHL2 = p(26);
Prod_of_ZEB = p(27);
Deg_of_ZEB = p(28);
Prod_of_miR200 = p(29);
Deg_of_miR200 = p(30);
Prod_of_SLUG = p(31);
Deg_of_SLUG = p(32);
GRHL2 = y(1);
ZEB = y(2);
miR200 = y(3);
SLUG = y(4);
HillsGRHL2_ZEB = Inh_of_GRHL2ToZEB*GRHL2^Num_of_GRHL2ToZEB/(GRHL2^Num_of_GRHL2ToZEB + Trd_of_GRHL2ToZEB^Num_of_GRHL2ToZEB)+ (1 -GRHL2^Num_of_GRHL2ToZEB/(GRHL2^Num_of_GRHL2ToZEB + Trd_of_GRHL2ToZEB^Num_of_GRHL2ToZEB));
HillsZEB_GRHL2 = Inh_of_ZEBToGRHL2*ZEB^Num_of_ZEBToGRHL2/(ZEB^Num_of_ZEBToGRHL2 + Trd_of_ZEBToGRHL2^Num_of_ZEBToGRHL2)+ (1 -ZEB^Num_of_ZEBToGRHL2/(ZEB^Num_of_ZEBToGRHL2 + Trd_of_ZEBToGRHL2^Num_of_ZEBToGRHL2));
HillsmiR200_ZEB = Inh_of_miR200ToZEB*miR200^Num_of_miR200ToZEB/(miR200^Num_of_miR200ToZEB + Trd_of_miR200ToZEB^Num_of_miR200ToZEB)+ (1 -miR200^Num_of_miR200ToZEB/(miR200^Num_of_miR200ToZEB + Trd_of_miR200ToZEB^Num_of_miR200ToZEB));
HillsZEB_miR200 = Inh_of_ZEBTomiR200*ZEB^Num_of_ZEBTomiR200/(ZEB^Num_of_ZEBTomiR200 + Trd_of_ZEBTomiR200^Num_of_ZEBTomiR200)+ (1 -ZEB^Num_of_ZEBTomiR200/(ZEB^Num_of_ZEBTomiR200 + Trd_of_ZEBTomiR200^Num_of_ZEBTomiR200));
HillsSLUG_miR200 = Inh_of_SLUGTomiR200*SLUG^Num_of_SLUGTomiR200/(SLUG^Num_of_SLUGTomiR200 + Trd_of_SLUGTomiR200^Num_of_SLUGTomiR200)+ (1 -SLUG^Num_of_SLUGTomiR200/(SLUG^Num_of_SLUGTomiR200 + Trd_of_SLUGTomiR200^Num_of_SLUGTomiR200));
HillsSLUG_ZEB = (Act_of_SLUGToZEB*SLUG^Num_of_SLUGToZEB/(SLUG^Num_of_SLUGToZEB + Trd_of_SLUGToZEB^Num_of_SLUGToZEB)+ (1 -SLUG^Num_of_SLUGToZEB/(SLUG^Num_of_SLUGToZEB + Trd_of_SLUGToZEB^Num_of_SLUGToZEB)))/Act_of_SLUGToZEB;
HillsZEB_ZEB = (Act_of_ZEBToZEB*ZEB^Num_of_ZEBToZEB/(ZEB^Num_of_ZEBToZEB + Trd_of_ZEBToZEB^Num_of_ZEBToZEB)+ (1 -ZEB^Num_of_ZEBToZEB/(ZEB^Num_of_ZEBToZEB + Trd_of_ZEBToZEB^Num_of_ZEBToZEB)))/Act_of_ZEBToZEB;
HillsSLUG_SLUG = (Act_of_SLUGToSLUG*SLUG^Num_of_SLUGToSLUG/(SLUG^Num_of_SLUGToSLUG + Trd_of_SLUGToSLUG^Num_of_SLUGToSLUG)+ (1 -SLUG^Num_of_SLUGToSLUG/(SLUG^Num_of_SLUGToSLUG + Trd_of_SLUGToSLUG^Num_of_SLUGToSLUG)))/Act_of_SLUGToSLUG;
dydt = [Prod_of_GRHL2*HillsZEB_GRHL2 - Deg_of_GRHL2*GRHL2
Prod_of_ZEB*HillsGRHL2_ZEB*HillsmiR200_ZEB*HillsSLUG_ZEB*HillsZEB_ZEB - Deg_of_ZEB*ZEB
Prod_of_miR200*HillsZEB_miR200*HillsSLUG_miR200 - Deg_of_miR200*miR200
Prod_of_SLUG*HillsSLUG_SLUG - Deg_of_SLUG*SLUG];
