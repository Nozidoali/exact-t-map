module top(a, b, c, y0);
input a, b, c;
output y0;
wire n0, n1, n2, n3;
assign n0 = a & b;
assign n1 = a & c;
assign n2 = b & c;
assign n3 = n0 | n1;
assign y0 = n3 | n2;
endmodule
