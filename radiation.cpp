double radiation(double tg, double n, double rad, float alpha, float kappa)
{
	return -(alpha*rad - kappa*(tg*tg))*n;
}
