function b = DCM2Quaternion(C)
% Converts DCM to quaternions via the Stanley Method

B2(1) = (1+trace(C))/4;
B2(2) = (1+2*C(1,1)-trace(C))/4;
B2(3) = (1+2*C(2,2)-trace(C))/4;
B2(4) = (1+2*C(3,3)-trace(C))/4;

[~,i] = max(B2);
switch i
	case 1
		b(1) = sqrt(B2(1));
		b(2) = (C(2,3)-C(3,2))/4/b(1);
		b(3) = (C(3,1)-C(1,3))/4/b(1);
		b(4) = (C(1,2)-C(2,1))/4/b(1);
	case 2
		b(2) = sqrt(B2(2));
		b(1) = (C(2,3)-C(3,2))/4/b(2);
		if (b(1)<0)
			b(2) = -b(2);
			b(1) = -b(1);
		end
		b(3) = (C(1,2)+C(2,1))/4/b(2);
		b(4) = (C(3,1)+C(1,3))/4/b(2);
	case 3
		b(3) = sqrt(B2(3));
		b(1) = (C(3,1)-C(1,3))/4/b(3);
		if (b(1)<0)
			b(3) = -b(3);
			b(1) = -b(1);
		end
		b(2) = (C(1,2)+C(2,1))/4/b(3);
		b(4) = (C(2,3)+C(3,2))/4/b(3);
	case 4
		b(4) = sqrt(B2(4));
		b(1) = (C(1,2)-C(2,1))/4/b(4);
		if (b(1)<0)
			b(4) = -b(4);
			b(1) = -b(1);
		end
		b(2) = (C(3,1)+C(1,3))/4/b(4);
		b(3) = (C(2,3)+C(3,2))/4/b(4);
end
b = b';
end