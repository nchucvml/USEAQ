#pragma once

#ifndef ARRAYLIST_H
#define ARRAYLIST_H

class ArrayList{
public:

	ArrayList(int cap)
	{
		this->cap = cap;
		A = new int[cap];
		max = 0;	
	};
	~ArrayList()
	{
		delete [] A;
	};


	void add(int item)
	{
		if(max<cap)
		{
			A[max] = item;
		}
		else
		{	
			cap = cap*2;
			int* ACopy = new int[cap];
			for(int i=0; i<max; i++)
			{
				ACopy[i] = A[i];
			}
			delete [] A;	
			A = ACopy;	
			A[max] = item;
		}
		max++;	
	};	

	void pop()
	{
		if(max > 1)
         max--;
	}

	ArrayList* clone()
	{
		ArrayList* clone = new ArrayList(cap);
		for(int i=0;i<max; i++)
			clone->add(A[i]);
		return clone;	
	};	

	double get(int index)
	{
		return A[index];	
	};
	void print()
	{
		for(int i=0; i<max; i++)
			std::cout << A[i] << " "; 
	};	
	void set (int index, int item)
	{	
		A[index] = item;
	};	
	int size()
	{
		return max;
	};
	bool operator==(ArrayList& other)
	{
		bool check = false;
		if(max != other.size()){
			return check;
		}else{
			for(int i=0; i<max; i++){
				if(other.get(i) != get(i))
				{
					return check;
				}
			}
			check = true;
			return check;
		}
	};

public:
	int cap,max;
	int* A;
};

#endif