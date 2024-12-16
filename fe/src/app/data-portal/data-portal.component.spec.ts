import { ComponentFixture, TestBed } from '@angular/core/testing';

import { DataPortalComponent } from './data-portal.component';

describe('DataPortalComponent', () => {
  let component: DataPortalComponent;
  let fixture: ComponentFixture<DataPortalComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [DataPortalComponent]
    })
    .compileComponents();

    fixture = TestBed.createComponent(DataPortalComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
