module top(x0, x1, x2, x3, x4, y0);
input x0, x1, x2, x3, x4;
output y0;
wire t0, t1, t2, t3, t4, t5, t6;
wire s0, s1, s2, s3, s4;
assign t0 = x0 & x1;
assign t1 = x0 & x3;
assign t2 = x1 & x2;
assign t3 = x1 & x3;
assign t4 = x1 & x4;
assign t5 = x2 & x3;
assign t6 = x3 & x4;
assign s0 = t0 ^ t1;
assign s1 = t2 ^ t3;
assign s2 = t4 ^ t5;
assign s3 = s0 ^ s1;
assign s4 = s2 ^ t6;
assign y0 = s3 ^ s4;
endmodule
