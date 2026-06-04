function x_next = boost_converter(x, u, T, L, C, R, Vin)
    % Modello non-lineare esatto al discreto del DC-DC Boost Converter
    x_next = zeros(2, 1);
    x_next(1) = x(1) + (T/L)*x(2) - (T/L)*(x(2) - Vin)*u;
    x_next(2) = -(T/C)*x(1) + (T/C)*x(1)*u + (1 - T/(R*C))*x(2);
end