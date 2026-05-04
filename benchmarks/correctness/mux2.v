module top(s, a, b, y0);
input s, a, b;
output y0;
wire n0, n1, n2;
assign n0 = s & a;
assign n1 = ~s;
assign n2 = n1 & b;
assign y0 = n0 | n2;
endmodule
