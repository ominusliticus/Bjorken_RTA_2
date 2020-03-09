#pragma once

#include <iostream>
#include <string>

using namespace std;

inline void fatalError(string message)
{
	cout << message << endl;
	char c;
	cout << "Terminating program. . ." << endl;
	cout << "Press any key to exit: " << endl;
	cin >> c;
	return;
}