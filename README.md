# CPP_REFPROP
修改NIST REFPROP自带的C++示例代码，方便批量调用refprop获取物性，例如密度，焓等

## 编译说明
将NIST REFPROP自带的refprop.dll 和 包含物质的两个文件夹fluids, mixtures和代码refprop.cpp, refprop.h放在同一目录下
## 代码说明
refprop中已经包含main函数，并有简单示例，参考修改即可，调用的函数返回指针，示例代码只获取了密度，也可获取其他参数

```
int main()
{
	Substance Nitrogen;
	char* fldname = "nitrogen.fld";
	Nitrogen.setSubstance(fldname);
	Nitrogen.showSubstance();
	double* pd  = Nitrogen.byTP(77,2000,1);
	printf(" density is: %f(mol/L).\n",pd[0]);

	for (int px = 1000; px < 2000; px = px + 100)
	{
		pd = Nitrogen.byTP(77.0, px*1.0, 1);
		printf("density is 77 (K)，pressure is %d (KPa)， density = %f6.2 (mol/L).\n",px,pd[0]);
	}
	return 0;
}
```
