OC = 0x00
FC = 0x10
# B_N = 0x01
# B_S = 0x02
# B_W = 0x04
# B_E = 0x08
# B_NW = 0x05
# B_SW = 0x06
# B_SE = 0x0A
# B_NE = 0x09
B_N = 0x01
B_S = 0x02
B_W = 0x04
B_E = 0x08
B_NW = B_N + B_W
B_SW = B_S + B_W
B_SE = B_S + B_E
B_NE = B_N + B_E

flag_dict = [("OC", OC),
         ("FC", FC),
         ("B_N", B_N),
         ("B_S", B_S),
         ("B_W", B_W),
         ("B_E", B_E),
         ("B_NW", B_NW),
         ("B_SW", B_SW),
         ("B_SE", B_SE),
         ("B_NE", B_NE)]
