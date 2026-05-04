module top(x0, x1, x2, x3, y0);
input x0, x1, x2, x3;
output y0;
wire t0, t1, t2, t3;
assign t0 = x0 & x1;
assign t1 = x0 & x3;
assign t2 = x1 & x2;
assign t3 = x2 & x3;
wire s1;
assign s1 = t0 ^ t1;
wire s2;
assign s2 = s1 ^ t2;
wire s3;
assign s3 = s2 ^ t3;
assign y0 = s3;
endmodule
